"""Pipeline progress: per-step status display and granular Snakemake resume targets."""
from __future__ import annotations

import json
import re
from pathlib import Path

from .plan_assembler import _stages_for_branch
from .provenance import current_branch
from .run_paths import RunPaths

# Ordered rows for status display (keys are monitor IDs, not always Snakemake stage IDs).
MONITOR_PIPELINE: tuple[str, ...] = (
    "plan_review",
    "s1a_ambient",
    "s1_rna_qc",
    "s2_atac_qc",
    "s3_doublets",
    "post_qc_review",
    "s4_rna_norm",
    "s5_atac_spectral",
    "s6_neighbors",
    "s7_sweep",
    "resolution_review",
    "s7_labels",
    "s8_umap",
)

RESOLUTION_REVIEW_STAGE = "s7_clustering"

HUMAN_GATES: frozenset[str] = frozenset({
    "plan_review", "post_qc_review", RESOLUTION_REVIEW_STAGE,
})

MONITOR_LABELS: dict[str, str] = {
    "plan_review": "plan_review",
    "s1a_ambient": "S1a",
    "s1_rna_qc": "S1",
    "s2_atac_qc": "S2",
    "s3_doublets": "S3",
    "post_qc_review": "qc_review",
    "s4_rna_norm": "S4",
    "s5_atac_spectral": "S5",
    "s6_neighbors": "S6",
    "s7_sweep": "S7",
    "resolution_review": "resolution_review",
    "s7_labels": "S7-labels",
    "s8_umap": "S8",
}

MONITOR_TASKS: dict[str, str] = {
    "plan_review": "Plan review",
    "s1a_ambient": "Ambient RNA correction",
    "s1_rna_qc": "RNA QC filtering",
    "s2_atac_qc": "ATAC QC filtering",
    "s3_doublets": "Doublet removal",
    "post_qc_review": "QC review",
    "s4_rna_norm": "RNA normalization",
    "s5_atac_spectral": "ATAC spectral embedding",
    "s6_neighbors": "PCA (RNA) + neighbor graph",
    "s7_sweep": "Resolution sweep",
    "resolution_review": "Resolution review",
    "s7_labels": "Cluster label assignment",
    "s8_umap": "UMAP and final outputs",
}

EXECUTE_MARKERS: dict[str, str] = {
    "s1a_ambient": "rna_decontaminated.h5ad",
    "s1_rna_qc": "qc_summary.json",   # persists after post_qc_review cleanup
    "s2_atac_qc": "qc_summary.json",  # persists after post_qc_review cleanup
    "s3_doublets": "calls.parquet",
    "s4_rna_norm": "rna_norm.h5ad",
    "s5_atac_spectral": "spectral_summary.json",
    "s6_neighbors": "rna_neighbors.h5ad",
    "s7_clustering": "rna_clustered.h5ad",
    "s8_umap": "s8_done.txt",
}

S7_SWEEP_MARKER = "sweep.parquet"

QC_EXECUTE_STAGES: tuple[str, ...] = (
    "s1a_ambient", "s1_rna_qc", "s2_atac_qc", "s3_doublets",
)
POST_QC_EXECUTE_STAGES: tuple[str, ...] = (
    "s4_rna_norm", "s5_atac_spectral", "s6_neighbors",
)

_RULE_ERROR_RE = re.compile(r"Error in rule ([A-Za-z0-9_]+):")

_FAILURE_MARKERS: tuple[str, ...] = (
    "RuleException:",
    "WorkflowError:",
    "Exiting because a job execution failed",
    "At least one job did not complete successfully",
)

_BLOCKABLE_STATES: frozenset[str] = frozenset({"pending", "in_progress", "cancelled"})

_KILL_ACTION_RE = re.compile(r"## Kill Action\s+```json\s+(.*?)\s+```", re.DOTALL)

# Snakemake rule name -> monitor row id (see MONITOR_PIPELINE).
_RULE_TO_MONITOR_ID: dict[str, str] = {
    "s1a_ambient_execute": "s1a_ambient",
    "s1_rna_qc_execute": "s1_rna_qc",
    "s2_atac_qc_execute": "s2_atac_qc",
    "s3_doublets_execute": "s3_doublets",
    "s4_rna_norm_execute": "s4_rna_norm",
    "s5_atac_spectral_execute": "s5_atac_spectral",
    "s6_neighbors_execute": "s6_neighbors",
    "s7_clustering_propose": "s7_sweep",
    "s7_clustering_execute": "s7_labels",
    "s8_umap_execute": "s8_umap",
    "plan_review_propose": "plan_review",
    "post_qc_review_propose": "post_qc_review",
}


def monitor_label(monitor_id: str) -> str:
    return MONITOR_LABELS.get(monitor_id, monitor_id)


def monitor_task(monitor_id: str) -> str:
    return MONITOR_TASKS.get(monitor_id, monitor_id)


def snakemake_rules_for_monitor(monitor_id: str) -> tuple[str, ...]:
    """Snakemake rule names whose logs indicate success/failure for a monitor row."""
    if monitor_id == "s7_sweep":
        return ("s7_clustering_propose",)
    if monitor_id == "s7_labels":
        return ("s7_clustering_execute",)
    if monitor_id == "resolution_review":
        return ()
    if monitor_id in HUMAN_GATES:
        return (f"{monitor_id}_propose",)
    if monitor_id in EXECUTE_MARKERS:
        return (f"{monitor_id}_propose", f"{monitor_id}_execute")
    return ()


def _branch_stages(paths: RunPaths) -> set[str]:
    return _stages_for_branch(current_branch(paths.parameters_yaml))


def _s7_applies(branch_stages: set[str]) -> bool:
    return RESOLUTION_REVIEW_STAGE in branch_stages


def _applies(monitor_id: str, branch_stages: set[str]) -> bool:
    if monitor_id in ("s7_sweep", "resolution_review", "s7_labels"):
        return _s7_applies(branch_stages)
    if monitor_id in HUMAN_GATES:
        return True
    return monitor_id in branch_stages


def execute_artifact(paths: RunPaths, stage: str) -> Path:
    return paths.artifact(stage, EXECUTE_MARKERS[stage])


def execute_done(paths: RunPaths, stage: str) -> bool:
    return execute_artifact(paths, stage).exists()


def _planning_done(paths: RunPaths) -> bool:
    """True when the merged planning job (s0_ingest_execute) has fully completed.

    The merged job loads + validates, assembles the plan, and runs the QC
    exploration in one shot. ``qc_explore.json`` is its last write, so requiring
    both the validation report and the explore JSON catches a job that died after
    ingest but before exploration. A dedicated predicate is used (not
    EXECUTE_MARKERS) because the explore marker lives in a different artifact dir
    than the stage id.
    """
    return (
        paths.artifact("s0_ingest", "validation_report.json").exists()
        and paths.artifact("qc_explore", "qc_explore.json").exists()
    )


def _s7_sweep_done(paths: RunPaths) -> bool:
    return paths.artifact(RESOLUTION_REVIEW_STAGE, S7_SWEEP_MARKER).exists()


def _log_indicates_failure(text: str) -> bool:
    return "Error in rule " in text or any(marker in text for marker in _FAILURE_MARKERS)


def _read_log_text(path: Path, *, max_bytes: int = 512_000) -> str:
    try:
        size = path.stat().st_size
    except OSError:
        return ""
    if size <= max_bytes:
        return path.read_text(errors="replace")
    with path.open("rb") as handle:
        handle.seek(max(0, size - max_bytes))
        return handle.read().decode(errors="replace")


def _snakemake_root(paths: RunPaths) -> Path:
    return paths.snakemake_workdir / ".snakemake"


def _cluster_logs_root(paths: RunPaths) -> Path:
    return _snakemake_root(paths) / "slurm_logs"


def _latest_rule_log(paths: RunPaths, rule_name: str) -> Path | None:
    log_dir = _cluster_logs_root(paths) / f"rule_{rule_name}"
    if not log_dir.is_dir():
        return None
    job_logs = list(log_dir.glob("*.log"))
    if not job_logs:
        return None
    return max(job_logs, key=lambda p: p.stat().st_mtime)


def _failed_rules_from_main_log(paths: RunPaths) -> set[str]:
    log_dir = _snakemake_root(paths) / "log"
    if not log_dir.is_dir():
        return set()
    logs = sorted(log_dir.glob("*.snakemake.log"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not logs:
        return set()
    text = _read_log_text(logs[0])
    if not _log_indicates_failure(text):
        return set()
    return set(_RULE_ERROR_RE.findall(text))


def _rule_is_failed(paths: RunPaths, rule_name: str, failed_rules: frozenset[str]) -> bool:
    if rule_name not in failed_rules:
        return False
    latest = _latest_rule_log(paths, rule_name)
    if latest is not None:
        return _log_indicates_failure(_read_log_text(latest))
    return True


def collect_failed_snakemake_rules(paths: RunPaths) -> frozenset[str]:
    """Snakemake rules whose logs report failure (per-rule cluster logs + main workflow log)."""
    failed: set[str] = set()

    cluster_root = _cluster_logs_root(paths)
    if cluster_root.is_dir():
        for rule_dir in cluster_root.iterdir():
            if not rule_dir.is_dir() or not rule_dir.name.startswith("rule_"):
                continue
            rule_name = rule_dir.name.removeprefix("rule_")
            latest = _latest_rule_log(paths, rule_name)
            if latest is not None and _log_indicates_failure(_read_log_text(latest)):
                failed.add(rule_name)

    failed.update(_failed_rules_from_main_log(paths))
    return frozenset(failed)


def _monitor_outputs_done(paths: RunPaths, monitor_id: str) -> bool:
    if monitor_id == "s7_sweep":
        return _s7_sweep_done(paths)
    if monitor_id == "s7_labels":
        return execute_done(paths, RESOLUTION_REVIEW_STAGE)
    if monitor_id in HUMAN_GATES or monitor_id == "resolution_review":
        return False
    return execute_done(paths, monitor_id)


def _child_job_to_rule(paths: RunPaths, job_id: str) -> str | None:
    """Map a SLURM child job id to the Snakemake rule name from cluster logs."""
    root = _cluster_logs_root(paths)
    if not root.is_dir():
        return None
    for rule_dir in root.iterdir():
        if not rule_dir.is_dir() or not rule_dir.name.startswith("rule_"):
            continue
        if (rule_dir / f"{job_id}.log").exists():
            return rule_dir.name.removeprefix("rule_")
    return None


def _load_monitor_kill_action(paths: RunPaths) -> dict | None:
    """Return the HPC monitor kill record written by Execution-MuAgent.

    Reads the structured ``kill_action`` field from latest_snapshot.json (the single
    machine contract). Falls back to parsing the ``## Kill Action`` block in
    latest_report.md only when the snapshot has no ``kill_action`` key — transitional
    support for a snapshot written by an older daemon.
    """
    snapshot = load_hpc_monitor_state(paths)
    if snapshot is not None and "kill_action" in snapshot:
        payload = snapshot.get("kill_action")
        if not isinstance(payload, dict) or not payload.get("attempted"):
            return None
        return payload

    # transitional: pre-structured-snapshot daemon — parse the debug markdown.
    report_path = paths.run_dir / "internal" / "hpc_monitor" / "latest_report.md"
    if not report_path.is_file():
        return None
    match = _KILL_ACTION_RE.search(report_path.read_text(errors="replace"))
    if not match:
        return None
    try:
        payload = json.loads(match.group(1))
    except json.JSONDecodeError:
        return None
    if not payload.get("attempted"):
        return None
    return payload


def collect_monitor_killed_monitor_ids(paths: RunPaths) -> frozenset[str]:
    """Monitor row ids whose SLURM child jobs were cancelled by the HPC monitor."""
    kill = _load_monitor_kill_action(paths)
    if kill is None:
        return frozenset()
    killed: set[str] = set()
    for action in kill.get("actions") or []:
        if action.get("role") != "child":
            continue
        job_id = str(action.get("job_id") or "").strip()
        if not job_id:
            continue
        rule = _child_job_to_rule(paths, job_id)
        if rule is None:
            continue
        monitor_id = _RULE_TO_MONITOR_ID.get(rule)
        if monitor_id is not None:
            killed.add(monitor_id)
    return frozenset(killed)


def _monitor_cancelled(
    paths: RunPaths,
    monitor_id: str,
    killed_monitor_ids: frozenset[str],
) -> bool:
    if monitor_id not in killed_monitor_ids:
        return False
    return not _monitor_outputs_done(paths, monitor_id)


def _monitor_failed(
    paths: RunPaths,
    monitor_id: str,
    failed_rules: frozenset[str],
) -> bool:
    if not failed_rules or _monitor_outputs_done(paths, monitor_id):
        return False
    for rule in snakemake_rules_for_monitor(monitor_id):
        if _rule_is_failed(paths, rule, failed_rules):
            return True
    return False


def _human_gate_state(
    paths: RunPaths,
    stage: str,
    failed_rules: frozenset[str],
    killed_monitor_ids: frozenset[str],
) -> str:
    if paths.approved_sentinel(stage).exists():
        return "approved"
    if paths.awaiting_sentinel(stage).exists():
        return "awaiting_approval"
    if stage == "post_qc_review" and paths.qc_review_summary_md.exists():
        return "awaiting_approval"
    if stage == RESOLUTION_REVIEW_STAGE and _s7_sweep_done(paths):
        return "awaiting_approval"
    if stage == "plan_review" and paths.plan_review_md.exists():
        return "awaiting_approval"
    if _monitor_failed(paths, stage, failed_rules):
        return "failed"
    if _monitor_cancelled(paths, stage, killed_monitor_ids):
        return "cancelled"
    return "pending"


def _automated_state(
    paths: RunPaths,
    monitor_id: str,
    stage_id: str,
    branch_stages: set[str],
    failed_rules: frozenset[str],
    killed_monitor_ids: frozenset[str],
) -> str:
    """Unified pending / in_progress / failed / cancelled / done / skipped for processing steps."""
    if stage_id not in branch_stages:
        return "skipped"
    if _monitor_outputs_done(paths, monitor_id):
        return "done"
    if _monitor_failed(paths, monitor_id, failed_rules):
        return "failed"
    if _monitor_cancelled(paths, monitor_id, killed_monitor_ids):
        return "cancelled"
    if paths.proposal(stage_id).exists():
        return "in_progress"
    return "pending"


def _monitor_state(
    paths: RunPaths,
    monitor_id: str,
    branch_stages: set[str],
    failed_rules: frozenset[str],
    killed_monitor_ids: frozenset[str],
) -> str:
    if monitor_id == "resolution_review":
        return _human_gate_state(paths, RESOLUTION_REVIEW_STAGE, failed_rules, killed_monitor_ids)

    if monitor_id == "s7_sweep":
        return _automated_state(
            paths, monitor_id, RESOLUTION_REVIEW_STAGE, branch_stages, failed_rules,
            killed_monitor_ids,
        )

    if monitor_id == "s7_labels":
        gate = _human_gate_state(paths, RESOLUTION_REVIEW_STAGE, failed_rules, killed_monitor_ids)
        if gate != "approved":
            return "pending"
        return _automated_state(
            paths, monitor_id, RESOLUTION_REVIEW_STAGE, branch_stages, failed_rules,
            killed_monitor_ids,
        )

    if monitor_id in HUMAN_GATES:
        return _human_gate_state(paths, monitor_id, failed_rules, killed_monitor_ids)

    return _automated_state(
        paths, monitor_id, monitor_id, branch_stages, failed_rules, killed_monitor_ids,
    )


def _apply_upstream_blocked(
    rows: list[tuple[str, str, str]],
) -> list[tuple[str, str, str]]:
    """After a failed or monitor-cancelled step, downstream pending rows become blocked."""
    blocked = False
    out: list[tuple[str, str, str]] = []
    for label, task, state in rows:
        if blocked and state in _BLOCKABLE_STATES:
            state = "blocked"
        if state in ("failed", "cancelled"):
            blocked = True
        out.append((label, task, state))
    return out


def stage_states(paths: RunPaths) -> list[tuple[str, str, str]]:
    """Return (short_label, task_name, state) for each applicable monitor row."""
    branch_stages = _branch_stages(paths)
    failed_rules = collect_failed_snakemake_rules(paths)
    killed_monitor_ids = collect_monitor_killed_monitor_ids(paths)
    rows = [
        (
            monitor_label(mid),
            monitor_task(mid),
            _monitor_state(paths, mid, branch_stages, failed_rules, killed_monitor_ids),
        )
        for mid in MONITOR_PIPELINE
        if _applies(mid, branch_stages)
    ]
    return _apply_upstream_blocked(rows)


def _last_incomplete_execute(
    paths: RunPaths,
    stages: tuple[str, ...],
    branch_stages: set[str],
) -> str | None:
    """Return the LAST incomplete stage in the phase as a Snakemake target.

    Targeting the last incomplete stage lets Snakemake build a single DAG that
    spans all remaining stages in the phase: it sees intermediate outputs are
    missing and chains s1a → s1 → s2 → s3 in one head-job submission.
    Targeting the *first* incomplete stage (the old behaviour) would submit only
    that single rule per head-job, requiring a separate submit call per stage.
    """
    last = None
    for stage in stages:
        if stage not in branch_stages:
            continue
        if not execute_done(paths, stage):
            last = f"{stage}_execute"
    return last


def infer_resume_target(run_dir: Path | str) -> str:
    """Pick the Snakemake target from the last incomplete step in the current phase.

    Using the last incomplete stage ensures Snakemake builds the full remaining
    phase DAG in a single head-job submission, rather than one stage per job.
    """
    paths = RunPaths(run_dir)
    branch_stages = _branch_stages(paths)

    # Planning phase: the merged s0_ingest job (load + validate + plan + QC
    # exploration) runs before checkpoint #1. Route it through `submit` so
    # Execution-MuAgent owns the cluster submission.
    if not _planning_done(paths):
        return "s0_ingest_execute"
    if not paths.approved_sentinel("plan_review").exists():
        return "plan_review_propose"

    incomplete = _last_incomplete_execute(paths, QC_EXECUTE_STAGES, branch_stages)
    if incomplete is not None:
        return incomplete
    if not paths.approved_sentinel("post_qc_review").exists():
        return "post_qc_review_propose"

    incomplete = _last_incomplete_execute(paths, POST_QC_EXECUTE_STAGES, branch_stages)
    if incomplete is not None:
        return incomplete
    if not paths.approved_sentinel(RESOLUTION_REVIEW_STAGE).exists():
        return "s7_clustering_propose"
    if not execute_done(paths, RESOLUTION_REVIEW_STAGE):
        return "s7_clustering_execute"
    if not execute_done(paths, "s8_umap"):
        return "s8_umap_execute"
    return "all"


def load_hpc_monitor_state(paths: RunPaths) -> dict | None:
    """Load latest_snapshot.json written by Execution-MuAgent.

    Returns the full JSON dict (including the nested "monitor_state" key when
    present) or None if the file is absent or unreadable.
    """
    snapshot_path = paths.run_dir / "internal" / "hpc_monitor" / "latest_snapshot.json"
    if not snapshot_path.is_file():
        return None
    try:
        return json.loads(snapshot_path.read_text(errors="replace"))
    except (json.JSONDecodeError, OSError):
        return None


def load_hpc_findings(paths: RunPaths) -> list[dict]:
    """Return the current-check findings list from latest_snapshot.json.

    Each entry is ``{"severity", "code", "message"}`` as written by the
    Execution-MuAgent daemon. Empty list when absent — Processing renders findings
    from this structured data, never from latest_report.md prose.
    """
    snapshot = load_hpc_monitor_state(paths)
    if not snapshot:
        return []
    findings = snapshot.get("findings")
    if not isinstance(findings, list):
        return []
    return [f for f in findings if isinstance(f, dict)]


def load_latest_hpc_submission(paths: RunPaths) -> dict | None:
    """Load latest_submission.json written by Execution-MuAgent."""
    path = paths.run_dir / "internal" / "hpc_monitor" / "latest_submission.json"
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(errors="replace"))
    except (json.JSONDecodeError, OSError):
        return None


def required_human_approvals(target: str) -> tuple[str, ...]:
    """Human checkpoint sentinels that must exist before running ``target``."""
    # Planning targets run BEFORE checkpoint #1 — they cannot require the
    # plan_review sentinel (which only exists after they produce the plan).
    # Without this short-circuit, `submit` would deadlock demanding an
    # impossible approval.
    if target in {"s0_ingest_execute", "plan_review_propose"}:
        return ()
    base = ("plan_review",)
    qc_targets = {f"{s}_execute" for s in QC_EXECUTE_STAGES} | {"post_qc_review_propose"}
    mid_targets = {f"{s}_execute" for s in POST_QC_EXECUTE_STAGES} | {"s7_clustering_propose"}
    if target in qc_targets:
        return base
    if target in mid_targets:
        return base + ("post_qc_review",)
    return base + ("post_qc_review", RESOLUTION_REVIEW_STAGE)
