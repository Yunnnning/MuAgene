"""User-facing plan review: concise summary built from P1 context + S0 ingest + P2 plan."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import provenance


def _load_json(p: Path) -> dict:
    return json.loads(p.read_text()) if p.exists() else {}


def build_summary(run_dir: Path | str) -> list[dict[str, Any]]:
    """Produce an ordered list of review items: {label, value, reason, certainty}."""
    from .run_paths import RunPaths
    paths = RunPaths(Path(run_dir))
    run_dir = paths.run_dir
    ctx = _load_json(paths.artifact("p1_context", "context_extraction.json"))
    ingest = _load_json(paths.artifact("s0_ingest", "validation_report.json"))
    plan = _load_json(paths.artifact("p2_plan", "preprocessing_plan.json"))

    def field(name: str) -> dict[str, Any]:
        return ctx.get("fields", {}).get(name, {})

    def param(stage: str, key: str) -> dict[str, Any]:
        return plan.get("stages", {}).get(stage, {}).get("parameters", {}).get(key, {})

    items: list[dict[str, Any]] = []

    # 1. Execution mode + HPC settings (from configure-execution / hpc.env)
    items.extend(_execution_review_items(run_dir, plan))

    # 2. Dataset type / modality
    modality = field("modality_type")
    items.append({
        "label": "Detected dataset type",
        "value": modality.get("value", "unknown"),
        "reason": modality.get("rationale", ""),
        "certainty": "certain" if modality.get("status") == "explicit" else "needs confirmation",
    })

    # 2. Organism + genome build
    organism = field("organism")
    genome = field("genome_build")
    items.append({
        "label": "Organism / genome build",
        "value": f"{organism.get('value', 'unknown')} / {genome.get('value', 'unknown')}",
        "reason": (f"{organism.get('rationale', '')} | {genome.get('rationale', '')}").strip(" |"),
        "certainty": "certain" if organism.get("status") == "explicit" and genome.get("confidence") != "low" else "needs confirmation",
    })

    # 3. Pairing
    pairing = ingest.get("pairing", {})
    items.append({
        "label": "Pairing (RNA ↔ ATAC)",
        "value": f"{pairing.get('status', 'unknown')} (overlap={pairing.get('overlap', 0):.2%})",
        "reason": f"{pairing.get('method', '')} at confidence={pairing.get('confidence', '')}; shared barcodes={pairing.get('n_shared', 0)}.",
        "certainty": "certain" if pairing.get("confidence") == "high" else "needs confirmation",
    })

    # 4. Ambient RNA correction (S1a) — confirm at plan review (10x: dataset + goal, not cells vs nuclei)
    ambient_method = param("s1a_ambient", "method")
    ambient_cap = param("s1a_ambient", "max_contamination")
    study_goal = param("s3_doublets", "study_goal")
    sample_type = field("sample_type").get("value", "unknown")
    rna_filtered_status = ingest.get("rna_filtered_status")
    has_raw = ingest.get("has_raw_matrix", False)
    if ambient_method.get("value") is not None:
        method_label = ambient_method.get("value", "auto")
        chosen_hint = ""
        if method_label == "auto":
            chosen_hint = (
                " → SoupX (raw + filtered both present)" if has_raw or rna_filtered_status == "raw"
                else " → DecontX (filtered only)"
            )
        elif method_label == "none":
            chosen_hint = " → pass-through (no correction)"
        goal_val = study_goal.get("value", "?")
        items.append({
            "label": "Ambient RNA correction (S1a)",
            "value": f"method={method_label}{chosen_hint}, cap={ambient_cap.get('value', '?')}",
            "reason": (
                f"study_goal={goal_val}, sample_type={sample_type}, "
                f"rna_filtered_status={rna_filtered_status}, has_raw={has_raw}. "
                "Default is auto unless run.yaml sets s1a_ambient_method. "
                "Confirm before S1a: use auto when rare populations or elevated ambient "
                "is likely; set method=none if you inspected Cell Ranger / markers and "
                "contamination is low (10x: not every dataset needs correction; nuclei "
                "can still have high ambient with debris). "
                f"Plan rationale: {ambient_method.get('rationale', '')}"
            ),
            "certainty": (
                "certain" if ambient_method.get("source") == "user" else "needs confirmation"
            ),
        })

    # 5. Key QC strategy
    sample_type = field("sample_type").get("value", "unknown")
    pct_mt_ceil = param("s1_rna_qc", "pct_mt_ceiling")
    pct_ribo_max = param("s1_rna_qc", "pct_ribo_max")
    tss_min = param("s2_atac_qc", "tss_enrichment_min")
    tss_max = param("s2_atac_qc", "tss_enrichment_max")
    nuc_max = param("s2_atac_qc", "nucleosome_signal_max")
    items.append({
        "label": "QC strategy",
        "value": (
            f"RNA: MAD on total_counts/n_genes, pct_mt ceiling={pct_mt_ceil.get('value', '?')}, "
            f"pct_ribo ceiling={pct_ribo_max.get('value', '?')} | "
            f"ATAC: TSS in ({tss_min.get('value', '?')}, {tss_max.get('value', '?')}), "
            f"MAD on log(n_fragments), nucleosome_signal<{nuc_max.get('value', '?')}"
        ),
        "reason": f"Sample type = {sample_type}; MAD thresholds adapt to the observed distribution.",
        "certainty": "certain",
    })

    # 5. Doublet policy — reconciliation is only meaningful when both detectors run.
    # Single-modality branches have one detector; the S3 item is informational, not
    # a checkpoint gate.
    branch = plan.get("workflow_branch", "paired")
    policy = param("s3_doublets", "removal_policy_recommendation")
    goal = param("s3_doublets", "study_goal")

    def _paired_doublet_policy_detail(recommended: str) -> list[str]:
        alt = "intersection" if recommended == "union" else "union"
        return [
            "- paired-multiome policy: Each modality runs its own doublet detector. "
            "You choose how to merge the two call lists before joint analysis:",
            "  - `union` — remove a cell if *either* modality flags it "
            "(stricter; more cells dropped).",
            "  - `intersection` — remove a cell only if *both* modalities flag it "
            "(more lenient; keeps cells where only one modality is suspicious).",
            f"- recommended: `{recommended}` — "
            + (
                "better when the priority is clean clustering and cell-type inference; "
                f"choose `{alt}` if retaining rare populations matters more than "
                "aggressively filtering ambiguous cells."
                if recommended == "union"
                else "better when retaining rare populations matters more than "
                "aggressively filtering ambiguous cells; "
                f"choose `{alt}` if the priority is clean clustering and cell-type inference."
            ),
        ]

    if branch == "paired":
        policy_val = policy.get("value", "?")
        items.append({
            "label": "Doublet removal policy",
            "value": policy_val,
            "detail": _paired_doublet_policy_detail(policy_val),
            "reason": "",
            "certainty": "needs confirmation",
        })
    elif branch == "separate":
        items.append({
            "label": "Doublet removal policy",
            "value": "independent (per-modality; Scrublet for RNA, SnapATAC2 for ATAC; raw calls preserved)",
            "reason": ("separate branch: modalities are independent samples with disjoint barcodes; "
                       "each modality's doublets removed by its own detector."),
            "certainty": "certain",
        })
    elif branch == "rna_only":
        items.append({
            "label": "Doublet removal policy",
            "value": "auto-applied (Scrublet only; raw calls preserved)",
            "reason": f"study_goal={goal.get('value', '?')}; single-detector branch — no reconciliation to confirm.",
            "certainty": "certain",
        })
    elif branch == "atac_only":
        items.append({
            "label": "Doublet removal policy",
            "value": "auto-applied (ATAC detector only; raw calls preserved)",
            "reason": f"study_goal={goal.get('value', '?')}; single-detector branch — no reconciliation to confirm.",
            "certainty": "certain",
        })

    # 6. Clustering
    grid = param("s7_clustering", "leiden_resolution_grid")
    rna_tilt = param("s7_clustering", "rna_tilt")
    atac_tilt = param("s7_clustering", "atac_tilt")
    items.append({
        "label": "Clustering strategy",
        "value": f"Leiden sweep {grid.get('value', [])} → stable-region knee (RNA={rna_tilt.get('value', '?')}, ATAC={atac_tilt.get('value', '?')})",
        "reason": "Per-modality resolution picked from a stability plateau, not a single-metric optimum.",
        "certainty": "certain",
    })

    # 7. Output location
    branch = plan.get("workflow_branch", "unknown")
    _s8 = paths.stage_dir("s8_umap")
    outputs = {
        "paired": f"{_s8}/processed.h5mu",
        "separate": f"{_s8}/rna_processed.h5ad + atac_processed.h5ad",
    }.get(branch, "(branch unknown)")
    items.append({
        "label": "Output location",
        "value": outputs,
        "reason": f"workflow_branch={branch}; run_manifest.json written alongside as handoff artifact.",
        "certainty": "certain",
    })

    # 8. Missing / uncertain
    missing = [k for k, v in ctx.get("fields", {}).items() if isinstance(v, dict) and v.get("status") == "missing"]
    low_conf = [k for k, v in ctx.get("fields", {}).items() if isinstance(v, dict) and v.get("confidence") == "low"]
    conflicts = ctx.get("conflicts", [])
    parts = []
    if missing:
        parts.append(f"missing: {', '.join(missing)}")
    if low_conf:
        parts.append(f"low-confidence: {', '.join(low_conf)}")
    if conflicts:
        parts.append(f"conflicts: {len(conflicts)}")
    items.append({
        "label": "Missing / uncertain info",
        "value": "; ".join(parts) if parts else "none",
        "reason": "Derived from P1 context field statuses and conflict list.",
        "certainty": "certain" if not parts else "needs confirmation",
    })

    return items


def _execution_review_items(run_dir: Path, plan: dict[str, Any]) -> list[dict[str, Any]]:
    from .run_paths import RunPaths

    exec_block = plan.get("execution") or {}
    mode = exec_block.get("mode", "local")
    settings = exec_block.get("settings") or {}
    params = provenance.load(RunPaths(run_dir).parameters_yaml)
    mode_entry = params.get("execution.mode") or {}
    user_configured = isinstance(mode_entry, dict) and mode_entry.get("source") == "user"

    items: list[dict[str, Any]] = [{
        "label": "Execution mode",
        "value": mode,
        "reason": exec_block.get(
            "s0_policy",
            "How heavy stages (S0 onward) are dispatched: local machine vs PBS/SLURM.",
        ),
        "certainty": "certain" if user_configured else "needs confirmation",
    }]

    if mode == "local":
        hpc_path = exec_block.get("hpc_env_path")
        items.append({
            "label": "HPC configuration",
            "value": hpc_path or "none (local default)",
            "reason": (
                "Cluster vars not required for a fully local run. "
                "If S0 fails locally, the agent may configure HPC and retry ingest on the cluster."
            ),
            "certainty": "certain" if not hpc_path else "needs confirmation",
        })
        return items

    # PBS or SLURM — surface scheduler settings for review
    parts: list[str] = []
    missing: list[str] = []
    if mode == "pbs":
        for label, key in (("queue", "pbs_queue"), ("project", "pbs_project")):
            val = settings.get(key)
            if val:
                parts.append(f"{label}={val}")
            elif key == "pbs_queue":
                missing.append("pbs_queue")
        for label, key in (("notify", "notify_email"), ("scale", "resources_scale"), ("conda", "conda_env")):
            val = settings.get(key)
            if val:
                parts.append(f"{label}={val}")
    else:
        for label, key in (("partition", "slurm_partition"), ("account", "slurm_account")):
            val = settings.get(key)
            if val:
                parts.append(f"{label}={val}")
            elif key == "slurm_partition":
                missing.append("slurm_partition")
        for label, key in (("notify", "notify_email"), ("scale", "resources_scale"), ("conda", "conda_env")):
            val = settings.get(key)
            if val:
                parts.append(f"{label}={val}")

    env_ref = exec_block.get("hpc_env_path") or "deliverables/pre_run/config/hpc.env"
    items.append({
        "label": "HPC configuration",
        "value": ", ".join(parts) if parts else "not set",
        "reason": (
            f"Source `{env_ref}` before cluster submit/resume. "
            + (f"Missing required: {', '.join(missing)}." if missing else "Scheduler settings recorded for this run.")
        ),
        "certainty": "needs confirmation" if missing else ("certain" if user_configured else "needs confirmation"),
    })
    return items


def _render_concise_section(items: list[dict[str, Any]]) -> str:
    lines: list[str] = []
    for it in items:
        tag = "✓" if it["certainty"] == "certain" else "?"
        lines.append(f"- **{it['label']}** [{tag} {it['certainty']}]")
        lines.append(f"  - value: `{it['value']}`")
        for detail_line in it.get("detail", []):
            lines.append(f"  {detail_line}")
        if it.get("reason"):
            lines.append(f"  - reason: {it['reason']}")
    return "\n".join(lines)


def render_summary_text(items: list[dict[str, Any]]) -> str:
    """Concise review bullets only (legacy helper for tests / partial display)."""
    return (
        "# Preprocessing plan review — summary\n\n"
        + _render_concise_section(items)
        + "\n"
    )


def render_merged_markdown(run_dir: Path | str) -> str:
    """Full plan-review deliverable: concise summary + parameter appendix."""
    from .plan_assembler import render_plan_appendix
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    paths = RunPaths(run_dir)
    items = build_summary(run_dir)
    plan = _load_json(paths.artifact("p2_plan", "preprocessing_plan.json"))
    parts = [
        "# Preprocessing plan review",
        "",
        "## Summary",
        "",
        _render_concise_section(items),
        "",
        "_Machine-readable plan: `internal/artifacts/p2_plan/preprocessing_plan.json` "
        "and `internal/parameters.yaml`. Full parameter listing in the appendix below._",
        "",
        "---",
        "",
    ]
    if plan:
        parts.append(render_plan_appendix(plan).rstrip())
        parts.append("")
    else:
        parts.append("_Appendix unavailable — `preprocessing_plan.json` not found._")
        parts.append("")
    return "\n".join(parts)


def write_summary(run_dir: Path | str) -> Path:
    """Write merged plan-review markdown and concise plan_summary.md."""
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    paths = RunPaths(run_dir)
    items = build_summary(run_dir)
    merged = render_merged_markdown(run_dir)
    paths.plan_review_md.parent.mkdir(parents=True, exist_ok=True)
    paths.plan_review_md.write_text(merged)
    paths.plan_summary_md.write_text(render_summary_text(items))
    return paths.plan_review_md
