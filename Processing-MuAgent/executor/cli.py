"""Processing-MuAgent CLI — wraps Snakemake for the interactive checkpointed workflow."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import click
import yaml

from . import approval, context as _ctx, hpc, plan_review as _pr, provenance, specs as _specs, stage_progress as _sp
from .log import log_event
from .run_paths import RunPaths


PACKAGE_DIR = Path(__file__).resolve().parent.parent  # Processing-MuAgent/
SNAKEFILE = PACKAGE_DIR / "workflow" / "Snakefile"

EXECUTOR_CHOICE = click.Choice(["local", "slurm"])
# Cluster-only executor for `submit` — `run` is local-only, so `local` is not a
# valid submit target (all cluster execution is owned by Execution-MuAgent).
CLUSTER_EXECUTOR_CHOICE = click.Choice(["slurm"])

# s0_ingest is the merged planning compute (load + validate + assemble plan +
# QC exploration); the former standalone p2_plan stage no longer exists.
STAGES = ["p1_context", "plan_review", "s0_ingest",
          "s1a_ambient", "s1_rna_qc", "s2_atac_qc", "s3_doublets",
          "post_qc_review",
          "s4_rna_norm", "s5_atac_spectral", "s6_neighbors", "s7_clustering", "s8_umap",
          "qc_handoff"]

HUMAN_CHECKPOINT_STAGES = ("plan_review", "post_qc_review")
STAGE_ALIASES = {
    "qc_review": "post_qc_review",
}
STAGE_DISPLAY = {
    "post_qc_review": "qc_review",
}
AUTOMATED_STAGES = tuple(s for s in STAGES if s not in HUMAN_CHECKPOINT_STAGES)


def _canonical_stage(stage: str) -> str:
    return STAGE_ALIASES.get(stage, stage)


def _display_stage(stage: str) -> str:
    return STAGE_DISPLAY.get(stage, stage)


def _resolve_run_dir(config_path: Path | str) -> Path:
    with Path(config_path).open() as f:
        cfg = yaml.safe_load(f) or {}
    rd = cfg.get("run_dir")
    if not rd:
        raise click.ClickException("run.yaml must set 'run_dir'")
    return Path(rd).expanduser().resolve()


@click.group()
def main() -> None:
    """Processing-MuAgent: multiome preprocessing subagent (stops after per-modality UMAP)."""


@main.command()
@click.option("--config", "config_path", required=True, type=click.Path())
def init(config_path: str) -> None:
    """Initialize a run directory.

    Creates the `internal/` and `deliverables/` scaffolds, copies the user's
    config into its canonical user-facing location `deliverables/plan/config/run.yaml`,
    and writes the Biological Context Report template into
    `deliverables/plan/config/biological_context.md`.
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    paths.ensure()
    # Config goes to its canonical deliverable location (Snakemake will read it
    # from there via --configfile — no separate internal copy).
    shutil.copy(config_path, paths.run_yaml)
    _ctx.write_template(paths.biological_context_md)
    click.echo(f"Initialized {run_dir}")
    click.echo(f"Fill {paths.biological_context_md} (optional but recommended).")


@main.command()
@click.argument("stage")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def propose(stage: str, config_path: str) -> None:
    """Run the <stage>_propose rule (local — propose rules are localrules)."""
    stage = _canonical_stage(stage)
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    _snakemake(["--configfile", str(paths.run_yaml), f"{stage}_propose"], run_dir)


# S1a recompute caches: regenerated whenever S1a re-runs, so they are both the
# stale-on-revise set (see _S1A_QC_ARTIFACTS) AND safe to delete once QC is
# approved (no further S1a re-run can occur). Single source of truth for both.
_S1A_REGEN_CACHES = [
    ("s1a_ambient", "tsne_coords_cache.parquet"),
    ("s1a_ambient", "cell_totals.parquet"),
]


def _run_config(run_dir: Path) -> dict:
    """Best-effort load of the run's canonical run.yaml (returns {} on any failure)."""
    try:
        return yaml.safe_load(RunPaths(run_dir).run_yaml.read_text()) or {}
    except Exception:
        return {}


def _retain_for_integration(run_dir: Path) -> bool:
    """Whether to KEEP the prepared ATAC fragment caches past the post_qc gate.

    Default True: Integration-MuAgent re-counts fragments against a consensus peak
    set, so the caches must survive approval. A single-sample-and-done run can set
    `retain_for_integration: false` in run.yaml to delete them and reclaim the disk.
    """
    return bool(_run_config(run_dir).get("retain_for_integration", True))


def _cleanup_qc_intermediates(run_dir: Path) -> list[str]:
    """Remove large QC-only working files after post_qc_review is approved.

    This is the single authority that deletes QC intermediate working files. None
    of these are declared Snakemake outputs — every QC/ingest stage declares only its
    durable marker (s0 validation_report.json, s1a summary.json, s1/s2 qc_summary.json)
    that carries the dependency edges, so removing them here does NOT make any rule's
    declared output "missing" and therefore never triggers a re-run of S0/S1a/S1/S2/S3
    on a later submit. None is read by any post-gate stage (S4–S8); the post-QC h5mu is
    the canonical store S4/S5 read. Removed:
      - rna_qc.h5ad / atac_qc.h5ad — untracked QC matrices, consumed only by S3.
      - atac_snap.h5ad — the pre-filter SnapATAC2 import.
      - atac_snap_explore.h5ad — the qc_explore ATAC import (reused by S2 during
        QC, no longer needed once QC is approved).
      - rna_ingest.h5ad / metadata_minimal.tsv — the S0 raw RNA ingest (~200 MB) and
        the unused reconstructed metadata TSV. rna_ingest.h5ad is consumed only by S1a
        (read by path); metadata_minimal.tsv has no reader. Dead once QC is approved.
      - rna_decontaminated.h5ad — the S1a ambient-corrected RNA (~400 MB), consumed
        only by S1 (read by path). Dead once QC is approved (S4 reads the post-QC h5mu).
      - atac_fragments_cbf[_chrnorm].tsv.gz (+ .tbi) — the chr-normalised fragment
        caches written by io.prepare_fragments_for_snapatac. These are the single
        biggest QC artifact and are RETAINED by default: Integration-MuAgent
        re-counts them against a consensus peak set (reading their *contents*, not
        just the cached filename recorded in parameters.yaml), so the old "dead once
        QC is approved" assumption no longer holds. They are deleted here only when
        run.yaml sets `retain_for_integration: false` (single-sample run reclaiming disk).
      - tsne_coords_cache.parquet / cell_totals.parquet — S1a recompute caches.

    Keeps the durable markers (validation_report.json — also read post-gate by S5 —,
    summary.json, qc_summary.json), the s1/s2 qc_metrics parquets (the durable
    QC-metrics record, consumed by the QC-review summary), and all S3+ artifacts so
    downstream stages are unaffected. Deleting an already-absent file is a no-op.

    Trade-off: deleting rna_ingest.h5ad / rna_decontaminated.h5ad means a *post-approval*
    QC `revise` of S1a/S1 can no longer re-run from these caches (they are gone and not
    a tracked input, so Snakemake will not regenerate them). Revising QC thresholds is a
    pre-approval activity; once the gate is approved, QC is committed.
    """
    rp = RunPaths(run_dir)
    # Untracked QC working matrices + the S0/S1a heavy RNA caches (DAG edges are the
    # durable markers validation_report.json / summary.json / qc_summary.json).
    targets = [
        rp.artifact("s1_rna_qc",  "rna_qc.h5ad"),
        rp.artifact("s2_atac_qc", "atac_qc.h5ad"),
        rp.artifact("s2_atac_qc", "atac_snap.h5ad"),
        rp.artifact("qc_explore", "atac_snap_explore.h5ad"),
        rp.artifact("s0_ingest",  "rna_ingest.h5ad"),
        rp.artifact("s0_ingest",  "metadata_minimal.tsv"),
        rp.artifact("s1a_ambient", "rna_decontaminated.h5ad"),
    ]
    # Chr-normalised fragment caches: both naming variants + tabix index, in
    # whichever stage dir they landed in (qc_explore import or an S2 re-derive).
    # RETAINED by default — Integration-MuAgent re-counts them against a consensus
    # peak set, so it reads their contents (not just the cached filename). Deleted
    # only when run.yaml opts out via `retain_for_integration: false`.
    if not _retain_for_integration(run_dir):
        for stage in ("qc_explore", "s2_atac_qc"):
            for name in ("atac_fragments_cbf_chrnorm.tsv.gz", "atac_fragments_cbf.tsv.gz"):
                targets.append(rp.artifact(stage, name))
                targets.append(rp.artifact(stage, name + ".tbi"))
    # S1a recompute caches (no S1a re-run can happen after approval).
    targets += [rp.artifact(s, f) for (s, f) in _S1A_REGEN_CACHES]

    deleted: list[str] = []
    for p in targets:
        if p.exists():
            p.unlink()
            deleted.append(str(p))
    return deleted


# S4–S8 working files that `finish-cleanup` deletes once the run's final processed
# deliverable exists. All are content-duplicates of the processed h5mu / h5ads
# (normalized RNA, PCA/neighbors, Leiden labels, ATAC spectral embedding, the
# exported peak/tile matrix + its peak-export scratch) — nothing downstream re-reads
# them. None is a declared Snakemake output, and every S4..S8 rule edge depends on the
# durable *_summary.json / s8_done.txt markers (see _PROCESS_MARKERS) instead, so
# deleting them never triggers a re-run on a later `submit --target all`.
_PROCESS_INTERMEDIATES = [
    # RNA chain (empty stubs on atac_only are removed too)
    ("s4_rna_norm", "rna_norm.h5ad"),
    ("s6_neighbors", "rna_neighbors.h5ad"),
    ("s7_clustering", "rna_clustered.h5ad"),
    # ATAC chain: S5 spectral working file + exported feature sidecars + peak-export
    # scratch (peak_matrix_*.h5ad / *_prepared.bed), and S7 ATAC labels.
    ("s5_atac_spectral", "atac_spectral.h5ad"),
    ("s5_atac_spectral", "feature_matrix.npz"),
    ("s5_atac_spectral", "feature_names.tsv"),
    ("s5_atac_spectral", "feature_kind.txt"),
    ("s5_atac_spectral", "peak_matrix_s2peaks.h5ad"),
    ("s5_atac_spectral", "peak_matrix_user.h5ad"),
    ("s5_atac_spectral", "_s2_peaks_prepared.bed"),
    ("s5_atac_spectral", "_user_peaks_prepared.bed"),
    ("s7_clustering", "atac_leiden_labels.parquet"),
]

# Durable per-stage done-markers that MUST survive finish-cleanup: `executor status`
# and the S4..S8 Snakemake edges key off them. (stage, filename).
_PROCESS_MARKERS = [
    ("s4_rna_norm", "norm_summary.json"),
    ("s5_atac_spectral", "spectral_summary.json"),
    ("s6_neighbors", "neighbors_summary.json"),
    ("s7_clustering", "clustering_summary.json"),
    ("s8_umap", "s8_done.txt"),
]


def _cleanup_process_intermediates(run_dir: Path) -> list[str]:
    """Remove the large S4–S8 working files once the processed deliverable exists.

    Mirrors `_cleanup_qc_intermediates` for the finish phase. Branch-awareness is by
    delete-if-exists: a branch simply never wrote the files it does not apply to
    (rna_only writes no ATAC sidecars), while an atac_only run's empty RNA stubs ARE
    removed. None of these is a declared Snakemake output — the durable markers in
    `_PROCESS_MARKERS` carry status + the DAG edges — so deletion never triggers a
    re-run. Deleting an absent file is a no-op. Returns the paths actually deleted.
    """
    rp = RunPaths(run_dir)
    deleted: list[str] = []
    for stage, name in _PROCESS_INTERMEDIATES:
        p = rp.artifact(stage, name)
        if p.exists():
            p.unlink()
            deleted.append(str(p))
    return deleted


def _processed_outputs_for_branch(rp: RunPaths, branch: str) -> list[Path]:
    """The final S8 processed deliverable(s) for `branch` — what finish-cleanup validates."""
    if branch == "paired":
        return [rp.processed_h5mu]
    out: list[Path] = []
    if branch in ("separate", "rna_only"):
        out.append(rp.rna_processed_h5ad)
    if branch in ("separate", "atac_only"):
        out.append(rp.atac_processed_h5ad)
    return out


def _s8_outputs_valid(run_dir: Path) -> tuple[bool, list[str]]:
    """Whether S8 produced its final processed deliverable(s) — the precondition for
    finish-cleanup. The processed h5mu / h5ads are the user-facing S8 output; a run
    that failed before producing them must keep its intermediates so it can resume.

    Primary check: the branch-derived processed file(s) exist and are non-empty.
    `run_manifest.json` is NOT required (a run can be S8-complete but not yet
    manifest-finalized — e.g. interrupted before the manifest rule); when present, the
    processed paths it records are additionally validated. Returns (ok, problems).
    """
    rp = RunPaths(run_dir)
    branch = provenance.current_branch(str(rp.parameters_yaml))
    problems: list[str] = []

    expected = _processed_outputs_for_branch(rp, branch)
    if not expected:
        problems.append(f"unknown workflow_branch {branch!r}; cannot locate S8 outputs")
    for p in expected:
        if not p.exists():
            problems.append(f"missing S8 output: {p}")
        elif p.stat().st_size == 0:
            problems.append(f"empty S8 output: {p}")

    manifest = rp.run_manifest_json
    if manifest.exists():
        import json
        try:
            outputs = (json.loads(manifest.read_text()) or {}).get("outputs") or {}
        except (json.JSONDecodeError, OSError) as e:
            problems.append(f"unreadable run_manifest.json: {e}")
            outputs = {}
        # outputs is a dict of {key: rel_path}; `figures` is a list (non-critical).
        for key, val in outputs.items():
            if key == "figures" or not isinstance(val, str):
                continue
            mp = run_dir / val
            if not mp.exists():
                problems.append(f"manifest output missing: {mp}")
            elif mp.is_file() and mp.stat().st_size == 0:
                problems.append(f"manifest output empty: {mp}")

    return (not problems, problems)


def _ensure_process_markers(run_dir: Path) -> list[str]:
    """Backfill any missing S4–S8 durable marker so `executor status` stays `done`
    after cleanup. Safe because the caller only invokes this once the processed
    deliverable is validated (the whole pipeline produced its final output), so every
    stage necessarily completed. Covers legacy / un-finalized runs predating the
    marker refactor (which lack norm/neighbors/clustering_summary.json and/or
    s8_done.txt). Returns the markers written.
    """
    import json
    rp = RunPaths(run_dir)
    written: list[str] = []
    for stage, name in _PROCESS_MARKERS:
        p = rp.artifact(stage, name)
        if p.exists():
            continue
        p.parent.mkdir(parents=True, exist_ok=True)
        if name.endswith(".json"):
            p.write_text(json.dumps({"stage": stage, "backfilled": True}, indent=2))
        else:
            p.write_text("backfilled\n")
        written.append(str(p))
    return written


# Per-stage QC artifacts that become stale when that stage is re-run, listed as
# (artifact_stage, filename). Revising a QC stage invalidates its own outputs and
# everything strictly downstream of it through S3 (the DAG is
# s1a_ambient -> s1_rna_qc -> s3_doublets and s2_atac_qc -> s3_doublets). The
# expensive chr-normalized fragment cache (atac_fragments_cbf_chrnorm.tsv.gz*) is
# intentionally NOT listed — it is reused across re-runs, and only deleted once QC
# is approved (by _cleanup_qc_intermediates, when no further re-run can occur).
_S1A_QC_ARTIFACTS = [
    ("s1a_ambient", "rna_decontaminated.h5ad"),  # untracked working file (read by S1 by path)
    ("s1a_ambient", "summary.json"),             # stage-done marker — delete to force re-run
] + _S1A_REGEN_CACHES
_S1_QC_ARTIFACTS = [
    ("s1_rna_qc", "rna_qc.h5ad"),
    ("s1_rna_qc", "qc_summary.json"),  # stage-done marker
]
_S2_QC_ARTIFACTS = [
    ("s2_atac_qc", "atac_qc.h5ad"),
    ("s2_atac_qc", "atac_snap.h5ad"),
    ("s2_atac_qc", "qc_summary.json"),  # stage-done marker
]
_S3_QC_ARTIFACTS = [
    ("s3_doublets", "rna_post_doublet.h5ad"),
    ("s3_doublets", "atac_post_doublet.h5ad"),
    ("s3_doublets", "calls.parquet"),
    ("s3_doublets", "joint_barcodes.txt"),
    ("s3_doublets", "overlap_summary.json"),
]
# stage -> stage-and-downstream artifact list (gate outputs always added on top).
_QC_INVALIDATION: dict[str, list[tuple[str, str]]] = {
    "s1a_ambient": _S1A_QC_ARTIFACTS + _S1_QC_ARTIFACTS + _S3_QC_ARTIFACTS,
    "s1_rna_qc":   _S1_QC_ARTIFACTS + _S3_QC_ARTIFACTS,
    "s2_atac_qc":  _S2_QC_ARTIFACTS + _S3_QC_ARTIFACTS,
    "s3_doublets": _S3_QC_ARTIFACTS,
}


def _qc_downstream_targets(run_dir: Path, stage: str) -> list[Path]:
    """The artifacts a revise of `stage` invalidates (whether or not they exist).

    The revised QC stage's own outputs, every strictly-downstream QC stage's
    outputs, AND the `post_qc_review_propose` gate outputs (proposal,
    awaiting_approval sentinel, qc_review_<run>.md, qc_summary_<run>.html). The
    gate outputs are always included: while they exist Snakemake reports "Nothing
    to be done" and silently skips the re-run even though upstream was invalidated.

    Pure — computes paths only, deletes nothing (so `revise --dry-run` can preview
    exactly what `revise` would remove).
    """
    if stage not in _QC_INVALIDATION:
        return []
    rp = RunPaths(run_dir)
    targets = [rp.artifact(s, f) for (s, f) in _QC_INVALIDATION[stage]]
    targets += [
        rp.proposal("post_qc_review"),
        rp.awaiting_sentinel("post_qc_review"),
        rp.qc_review_summary_md,
        rp.qc_summary_html,
        rp.post_qc_h5mu,
        rp.post_qc_manifest_json,
    ]
    return targets


def _invalidate_qc_downstream(run_dir: Path, stage: str) -> list[str]:
    """Delete the stale artifacts from `_qc_downstream_targets` so a re-run from
    `stage` actually re-executes. Deleting non-existent files is a no-op (safe to
    call before QC has ever run, e.g. a marker-gene revise at plan review).
    Returns the paths actually deleted (for transparent logging).
    """
    deleted: list[str] = []
    for p in _qc_downstream_targets(run_dir, stage):
        if p.exists():
            p.unlink()
            deleted.append(str(p))
    return deleted


def _echo_binding_constraint(paths: RunPaths, stage: str) -> None:
    """Print the current QC thresholds for a QC stage so the user can see which
    bound (MAD-derived vs floor vs ceiling) is actually binding before revising —
    the pct_mt floor-vs-ceiling case is the classic gotcha."""
    if stage not in ("s1_rna_qc", "s2_atac_qc"):
        return
    import json as _json
    qexp = paths.artifact("qc_explore", "qc_explore.json")
    if not qexp.exists():
        return
    try:
        th = (_json.loads(qexp.read_text()).get(stage) or {}).get("thresholds")
    except Exception:
        return
    if th:
        click.echo(f"  binding-constraint check ({stage}) — current thresholds "
                   "(compare MAD-derived vs floor/ceiling to see which is active):")
        for k, v in th.items():
            click.echo(f"    {k} = {v}")


def _revise_dry_run(run_dir: Path, paths: RunPaths, stage: str, key: str, value_parsed) -> None:
    """Preview a `revise`: the parameter change, the binding-constraint context, and
    the EXACT artifacts that would be deleted — mutating nothing. Closes the
    destructive-revise gap (revise at post_qc_review used to delete downstream
    outputs with no preview or undo)."""
    click.echo("DRY RUN — no changes made.")
    current = provenance.get_value(str(paths.parameters_yaml), key, None)
    cur_disp = repr(current) if current is not None else "(plan default; no override set)"
    click.echo(f"  param: {key}: {cur_disp} -> {value_parsed!r}")
    if not approval.is_approved(run_dir, "plan_review"):
        click.echo("  checkpoint: plan_review (not yet approved) -> would re-render the plan "
                   "deliverables (overlay); NO artifacts deleted.")
    else:
        existing = [p for p in _qc_downstream_targets(run_dir, stage) if p.exists()]
        if existing:
            click.echo(f"  checkpoint: post_qc_review -> WOULD DELETE {len(existing)} artifact(s):")
            for p in existing:
                click.echo(f"    - {p}")
            click.echo("  Re-running the affected stages regenerates them. Confirm before "
                       "running the real `revise`.")
        else:
            click.echo("  checkpoint: post_qc_review -> no existing downstream artifacts to delete.")
    _echo_binding_constraint(paths, stage)


def _regenerate_plan_deliverables(run_dir: Path) -> list[str]:
    """Re-render the plan-review deliverables after a `revise` at the (still
    unapproved) plan_review gate, so the overlay reflects the change.

    Cheap: the QC-exploration preview re-derives from the persisted per-cell
    metrics parquets (no h5ad reload, no fragment re-import), then the markdown
    and HTML re-render with the effective (override-overlaid) plan. No-op-safe
    when the metrics or plan are absent. Returns what was regenerated.
    """
    from executor import plan_review as _plan_review
    from executor import qc_explore

    paths = RunPaths(run_dir)
    regenerated: list[str] = []
    try:
        qc_explore.rederive_from_metrics(run_dir)
        regenerated.append("qc_explore preview")
    except Exception as exc:  # preview is best-effort; rendering proceeds either way
        log_event(run_dir, {"stage": "plan_review", "event": "qc_rederive_failed",
                            "error": str(exc)})
    try:
        _plan_review.write_summary(run_dir)
        _plan_review.write_plan_summary_html(run_dir)
        regenerated += [paths.plan_review_md.name, paths.plan_summary_html.name]
    except Exception as exc:
        log_event(run_dir, {"stage": "plan_review", "event": "plan_rerender_failed",
                            "error": str(exc)})
    return regenerated


_MARKER_GENE_GATE_MSG = (
    "Ambient RNA correction is planned but no marker genes are set and no explicit "
    "decision was recorded. Ask the user whether to check marker-gene expression "
    "before vs after correction (recommended, especially at elevated contamination) "
    "and to provide 5-10 gene symbols. Then re-run approve with one of:\n"
    "  - provide genes:  executor revise s1a_ambient s1a_ambient.marker_genes=\"[GENE1, GENE2]\" --config <cfg>\n"
    "  - defer to QC review:  approve plan_review --defer-marker-genes\n"
    "  - decline:  approve plan_review --skip-marker-genes\n"
    "Never invent or suggest gene symbols yourself."
)


def _apply_marker_gene_ack(run_dir: Path, ack: str | None) -> None:
    """Record an unattended-batch marker-gene decision (`--marker-genes defer|skip`)."""
    if not ack:
        return
    decision = {"defer": "deferred_to_qc", "skip": "declined"}[ack]
    _pr.record_marker_gene_decision(run_dir, decision)
    log_event(run_dir, {"stage": "plan_review", "event": "marker_gene_decision",
                        "decision": decision})


def _resolve_marker_gene_gate(run_dir: Path, *, defer: bool, skip: bool) -> None:
    """Enforce an explicit marker-gene decision before plan_review is approved."""
    if defer and skip:
        raise click.ClickException(
            "Pass at most one of --defer-marker-genes / --skip-marker-genes.")
    if defer or skip:
        decision = "deferred_to_qc" if defer else "declined"
        _pr.record_marker_gene_decision(run_dir, decision)
        log_event(run_dir, {"stage": "plan_review", "event": "marker_gene_decision",
                            "decision": decision})
        return
    if _pr.marker_gene_decision_pending(run_dir):
        raise click.ClickException(_MARKER_GENE_GATE_MSG)


@main.command()
@click.argument("stage")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--note", default="")
@click.option("--defer-marker-genes", "defer_marker_genes", is_flag=True,
              help="plan_review only: record an explicit choice to check marker "
                   "genes at QC review instead of now.")
@click.option("--skip-marker-genes", "skip_marker_genes", is_flag=True,
              help="plan_review only: record an explicit choice to decline the "
                   "before/after-ambient marker gene expression check.")
def approve(stage: str, config_path: str, note: str,
            defer_marker_genes: bool, skip_marker_genes: bool) -> None:
    """Write internal/checkpoints/<stage>.approved to unblock <stage>_execute."""
    stage = _canonical_stage(stage)
    run_dir = _resolve_run_dir(config_path)
    if stage == "plan_review":
        _resolve_marker_gene_gate(
            run_dir, defer=defer_marker_genes, skip=skip_marker_genes)
    elif defer_marker_genes or skip_marker_genes:
        raise click.ClickException(
            "--defer-marker-genes / --skip-marker-genes apply only to "
            "`approve plan_review`.")
    approval.approve(run_dir, stage, note=note)
    log_event(run_dir, {"stage": stage, "event": "approved", "note": note})
    if stage == "post_qc_review":
        deleted = _cleanup_qc_intermediates(run_dir)
        if deleted:
            log_event(run_dir, {"stage": "post_qc_review", "event": "qc_cleanup",
                                 "deleted": deleted})
            click.echo(f"Cleaned up {len(deleted)} intermediate QC object(s).")
    click.echo(f"Approved {_display_stage(stage)}")


@main.command(name="finish-cleanup")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def finish_cleanup(config_path: str) -> None:
    """Delete the large S4–S8 intermediate working files after the run is complete.

    Run this once S8 has produced the final processed deliverable. It first VALIDATES
    that the S8 output exists and is non-empty; if not, it refuses and leaves every
    intermediate in place so the pipeline can still resume from an intermediate stage.
    On success it backfills any missing durable stage markers (so `status` keeps
    reporting S4–S8 done) and removes the working h5ads/sidecars — content-duplicates
    of the processed deliverable. The deletions are not declared Snakemake outputs, so
    a later `submit --target all` does not re-run S4–S8.
    """
    run_dir = _resolve_run_dir(config_path)
    ok, problems = _s8_outputs_valid(run_dir)
    if not ok:
        raise click.ClickException(
            "S8 output not present/valid — refusing finish-cleanup so the run can "
            "resume from an intermediate step:\n  " + "\n  ".join(problems))
    backfilled = _ensure_process_markers(run_dir)
    if backfilled:
        log_event(run_dir, {"stage": "finish_cleanup", "event": "markers_backfilled",
                            "written": backfilled})
    deleted = _cleanup_process_intermediates(run_dir)
    if deleted:
        log_event(run_dir, {"stage": "finish_cleanup", "event": "process_cleanup",
                            "deleted": deleted})
        click.echo(f"Cleaned up {len(deleted)} S4–S8 intermediate object(s).")
    else:
        click.echo("No S4–S8 intermediates to clean.")


@main.command(name="qc-cleanup")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def qc_cleanup(config_path: str) -> None:
    """Delete the large QC/ingest intermediate working files of an approved run.

    This is the same cleanup `approve post_qc_review` runs automatically; expose it
    standalone so disk can be reclaimed on a run that was approved earlier (e.g. to
    apply an expanded cleanup set retroactively). REQUIRES `post_qc_review` to be
    approved — refuses otherwise, since the deleted caches (rna_ingest.h5ad,
    rna_decontaminated.h5ad, the QC matrices) are still needed while QC is in review.
    The durable markers survive, so nothing re-runs; deliverables are untouched.
    """
    run_dir = _resolve_run_dir(config_path)
    if not approval.is_approved(run_dir, "post_qc_review"):
        raise click.ClickException(
            "post_qc_review is not approved — refusing qc-cleanup. These caches are "
            "needed while QC is in review (and the before/after-ambient marker-gene "
            "check reads rna_decontaminated.h5ad). Approve QC first.")
    deleted = _cleanup_qc_intermediates(run_dir)
    if deleted:
        log_event(run_dir, {"stage": "post_qc_review", "event": "qc_cleanup",
                            "deleted": deleted})
        click.echo(f"Cleaned up {len(deleted)} intermediate QC object(s).")
    else:
        click.echo("No QC intermediates to clean.")


@main.command(name="declare-branch")
@click.argument("branch", type=click.Choice(["paired", "separate", "rna_only", "atac_only"]))
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def declare_branch(branch: str, config_path: str) -> None:
    """Declare the workflow branch up front (user assertion).

    Writes `plan.workflow_branch_declared` to parameters.yaml with source=user.
    S0 will confirm this matches its own detection, or raise with a clear diff.
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    paths.ensure()
    provenance.set_param(
        str(paths.parameters_yaml),
        "plan.workflow_branch_declared", branch,
        source="user", confidence="high",
        rationale=f"Declared via `executor declare-branch {branch}`.",
    )
    log_event(run_dir, {"stage": "declare_branch", "event": "declared", "branch": branch})
    click.echo(f"Declared workflow_branch={branch!r}; S0 will confirm at ingest time.")


@main.command(name="hpc-info")
def hpc_info() -> None:
    """Probe the login node for scheduler queues/partitions and current PMA_* env."""
    import json
    info = hpc.discover_site()
    click.echo(json.dumps(info, indent=2, sort_keys=True))


@main.command(name="configure-execution")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--mode", "mode", required=True, type=EXECUTOR_CHOICE,
              help="Execution backend: local or slurm.")
@click.option("--confirmed-by-user/--not-confirmed", "confirmed_by_user",
              default=False,
              help="Record that the USER explicitly confirmed this execution mode "
                   "(local vs HPC). `run`/`submit` refuse to launch any compute job "
                   "until this is set. Never pass it on the user's behalf without "
                   "having actually confirmed.")
@click.option("--slurm-partition", default=None, help="SLURM partition (PMA_SLURM_PARTITION).")
@click.option("--slurm-account", default=None, help="SLURM account (PMA_SLURM_ACCOUNT).")
@click.option("--resources-scale", default=None, type=float,
              help="Memory/walltime scale factor (PMA_RESOURCES_SCALE).")
@click.option("--conda-env", default=None, help="Conda env name for cluster jobs (PMA_CONDA_ENV).")
@click.option("--device", type=click.Choice(["cpu", "gpu"]), default="cpu", show_default=True,
              help="Compute device for GPU-capable stages (compute.device). 'gpu' is cluster-only "
                   "(--mode slurm + submit): routes those stages to the GPU partition/gres and "
                   "container; local mode (--mode local) is always CPU.")
@click.option("--gpu-partition", default=None,
              help="SLURM GPU partition for GPU-capable stages (PMA_SLURM_GPU_PARTITION).")
@click.option("--gpu-gres", default=None,
              help="SLURM GPU gres request, e.g. gpu:A5000:1 (PMA_SLURM_GPU_GRES). Run hpc-info to discover.")
@click.option("--gpu-conda-env", default=None,
              help="Conda env activated for GPU child jobs (PMA_CONDA_ENV_GPU), e.g. muagene-gpu.")
@click.option("--gpu-image", default=None,
              help="Machine-local path the GPU .sif is pulled to (default ~/.muagene/images/muagene-gpu.sif).")
@click.option("--gpu-image-uri", default=None,
              help="Pinned registry reference the GPU image is PULLED from, e.g. "
                   "docker://<registry>/muagene-gpu:<tag>. No machine builds the image locally. "
                   "Defaults to machine.config gpu_image_uri.")
@click.option("--singularity-module", default=None,
              help="Module to `module load` before singularity exec (e.g. singularityce/3.11.3). "
                   "Defaults to machine.config singularity_module.")
@click.option("--scratch", default=None,
              help="Optional node-local/fast scratch path to bind into the GPU container "
                   "(exported as PMA_GPU_BIND). The run directory and repo root are always "
                   "bound; this is for extra paths a stage writes outside the run dir.")
@click.option("--env-policy", type=click.Choice(["auto", "manual"]), default="auto", show_default=True,
              help="On a missing/stale env at submit: auto-provision (default) or fail loud with the command.")
def configure_execution(
    config_path: str,
    mode: str,
    confirmed_by_user: bool,
    slurm_partition: str | None,
    slurm_account: str | None,
    resources_scale: float | None,
    conda_env: str | None,
    device: str,
    gpu_partition: str | None,
    gpu_gres: str | None,
    gpu_conda_env: str | None,
    gpu_image: str | None,
    gpu_image_uri: str | None,
    singularity_module: str | None,
    scratch: str | None,
    env_policy: str,
) -> None:
    """Record execution mode and write deliverables/plan/config/site.config + hpc.env."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    paths.ensure()

    params_path = str(paths.parameters_yaml)
    prior_mode = provenance.get_value(params_path, "execution.mode", None)
    prior_confirmed = provenance.get_value(params_path, "execution.user_confirmed", False)

    provenance.set_param(
        params_path,
        "execution.mode", mode,
        source="user", confidence="high",
        rationale=f"Execution backend set via configure-execution --mode {mode}.",
    )

    if mode == "local" and device == "gpu":
        raise click.ClickException(
            "GPU is cluster-only: --device gpu requires --mode slurm (use submit). "
            "Local runs use --mode local with the default --device cpu.")

    # On HPC, gpu routes GPU-capable stages to the GPU partition/gres and container
    # (see workflow/resources.smk _GPU_CAPABLE). Stages that are not GPU-capable
    # always run on CPU regardless of this setting. Local mode is CPU-only.
    provenance.set_param(
        params_path,
        "compute.device", device,
        source="user", confidence="high",
        rationale=f"Compute device set via configure-execution --device {device}.",
    )

    # Explicit, auditable record of whether the USER confirmed this execution mode.
    # `run`/`submit` refuse to launch any compute job until this is true (see
    # `_enforce_execution_mode_gate`). Recording the mode alone is not enough —
    # the agent must never silently choose local vs HPC.
    #
    # Re-config semantics: an explicit --confirmed-by-user always confirms. Without
    # it, confirmation is PRESERVED only when the mode is unchanged (e.g. bumping
    # --resources-scale on the same backend) — so resource tweaks don't silently
    # un-confirm a run. A *changed* mode (or one never confirmed) resets to
    # unconfirmed, forcing a fresh user confirmation.
    if confirmed_by_user:
        confirmed = True
        confirmed_rationale = ("User explicitly confirmed the execution mode via "
                               "configure-execution --confirmed-by-user.")
    elif prior_confirmed and prior_mode == mode:
        confirmed = True
        confirmed_rationale = (f"Confirmation preserved across re-config of unchanged "
                               f"mode {mode!r} (no --confirmed-by-user needed for "
                               "resource-only changes).")
    else:
        confirmed = False
        confirmed_rationale = ("Execution mode recorded WITHOUT explicit user "
                               "confirmation (--not-confirmed / default, or mode "
                               "changed); run/submit will refuse compute.")
    provenance.set_param(
        params_path,
        "execution.user_confirmed", confirmed,
        source="user", confidence="high",
        rationale=confirmed_rationale,
    )
    if not confirmed:
        click.echo(
            "NOTE: execution mode recorded but NOT user-confirmed. `run`/`submit` "
            "will refuse to launch any compute job until you confirm local vs HPC "
            "with the user and re-run:\n"
            f"  Processing-MuAgent configure-execution --config {paths.run_yaml} "
            f"--mode {mode} --confirmed-by-user",
            err=True,
        )

    # Machine-level infra knobs default from ~/.muagene/machine.config (written once by
    # Execution-MuAgent `init-machine`), so the operator doesn't re-type manager/module/
    # image/env names per run. Precedence: explicit flag > machine.config > env var.
    mc = hpc.load_machine_config()
    settings: dict[str, str | None] = {
        "slurm_partition": slurm_partition or os.environ.get("PMA_SLURM_PARTITION"),
        "slurm_account": slurm_account or os.environ.get("PMA_SLURM_ACCOUNT"),
        "resources_scale": (
            str(int(resources_scale)) if resources_scale is not None
            else os.environ.get("PMA_RESOURCES_SCALE")
        ),
        "conda_env": (conda_env or mc.get("conda_env") or os.environ.get("PMA_CONDA_ENV")
                      or os.environ.get("CONDA_DEFAULT_ENV")),
        "device": device,
        "slurm_gpu_partition": gpu_partition or os.environ.get("PMA_SLURM_GPU_PARTITION"),
        "slurm_gpu_gres": gpu_gres or os.environ.get("PMA_SLURM_GPU_GRES"),
        "gpu_conda_env": gpu_conda_env or mc.get("gpu_conda_env") or os.environ.get("PMA_CONDA_ENV_GPU"),
        "gpu_image": gpu_image or mc.get("gpu_image") or os.environ.get("PMA_GPU_IMAGE"),
        "gpu_image_uri": gpu_image_uri or mc.get("gpu_image_uri") or os.environ.get("PMA_GPU_IMAGE_URI"),
        "singularity_module": (singularity_module or mc.get("singularity_module")
                               or os.environ.get("PMA_SINGULARITY_MODULE")),
        # Optional extra GPU-container bind (-> PMA_GPU_BIND). Flag > machine.config > env.
        "scratch": scratch or mc.get("scratch") or os.environ.get("PMA_GPU_BIND"),
        # Detected infra: None lets Execution auto-detect; machine.config pins them.
        "env_manager": mc.get("manager"),
        "container_runtime": mc.get("container_runtime"),
        "env_policy": env_policy or mc.get("policy") or "auto",
    }

    if mode == "local":
        click.echo(f"Execution mode: local (device={device}; no hpc.env written).")
        return

    if mode == "slurm" and not settings["slurm_partition"]:
        raise click.ClickException(
            "SLURM mode requires --slurm-partition or PMA_SLURM_PARTITION in the environment.")

    # GPU routing prerequisites (cluster-only; preprocessing stages are CPU-only —
    # _GPU_CAPABLE is empty until the integration subagent adds stages). Fail loud
    # rather than silently submitting with a misconfigured partition/env.
    if device == "gpu":
        click.echo(
            "NOTE: --device gpu prepares cluster GPU routing for the integration "
            "subagent (future). Processing-MuAgent preprocessing is CPU-only.",
            err=True,
        )
        # SLURM GPU is container-only: the job runs inside the PULLED image, so fail
        # loud now if the pinned image reference is missing rather than writing
        # image_uri: null and only discovering it at provision/submit
        # (gpu_image_unavailable).
        if mode == "slurm" and not settings["gpu_image_uri"]:
            raise click.ClickException(
                "SLURM --device gpu requires --gpu-image-uri — a pinned registry reference the GPU "
                "image is PULLED from (e.g. docker://<registry>/muagene-gpu:<tag>) — or gpu_image_uri "
                "in ~/.muagene/machine.config (set once via `Execution-MuAgent init-machine`). "
                "No machine builds the image locally.")
        if mode == "slurm" and not settings["slurm_gpu_gres"]:
            raise click.ClickException(
                "SLURM --device gpu requires --gpu-gres (e.g. gpu:A5000:1). Run `hpc-info` to discover "
                "the GPU partition/gres on this cluster.")

    site_cfg = hpc.write_site_config(paths.site_config, mode=mode, settings=settings)
    out = hpc.write_hpc_env(paths.hpc_env_sh, paths.site_config)
    log_event(run_dir, {"stage": "configure_execution", "event": "configured",
                        "mode": mode, "hpc_env": str(out), "site_config": str(site_cfg)})
    click.echo(f"Execution mode: {mode}")
    click.echo(f"Wrote {site_cfg}")
    click.echo(f"Wrote {out}  (derived from site.config)")
    click.echo("Source this file in your shell before submit/run on the cluster:")
    click.echo(f"  source {out}")


@main.command(name="regenerate-locks")
@click.option("--platform", "platforms", multiple=True, default=("linux-64",), show_default=True,
              help="conda platform(s) to lock for. MuAgene is linux-only; default linux-64.")
def regenerate_locks(platforms: tuple[str, ...]) -> None:
    """Regenerate the CPU conda-lock lockfile from workflow/envs/processing.yaml.

    The YAML is the human source of truth; the committed lock is what actually gets
    installed (solve-free, reproducible). Run this AFTER editing processing.yaml, then
    COMMIT the refreshed lock — `validate-env`/`submit` fail loud (`lock_stale_vs_yaml`)
    when the YAML's content hash no longer matches the lock's recorded `# source-sha256:`.
    Lock generation is a science-authoring act, so it lives in Processing-MuAgent.
    Requires conda-lock: `pip install 'Processing-MuAgent[dev]'`.
    """
    import hashlib
    import shutil
    import subprocess

    import yaml

    man = hpc.load_env_manifest()
    cpu = man.get("cpu") or {}
    yaml_path = hpc.REPO_ROOT / cpu["definition"]
    work = (hpc.REPO_ROOT / cpu["lock"]).parent
    if not yaml_path.exists():
        raise click.ClickException(f"CPU env YAML not found: {yaml_path}")
    # The CPU env is rendered with `conda-lock --kind explicit` — a conda-ONLY format that
    # silently drops any `pip:` subsection. A pip dep here would therefore never reach the
    # lock (so a freshly provisioned env would be missing it and fail validate-env at run
    # time, far from this command). Fail loud instead: every dependency must be a conda
    # package. `- pip` itself (a bare string, for `init-machine`'s editable agent installs)
    # is fine; only a `pip:` mapping is rejected.
    spec = yaml.safe_load(yaml_path.read_text()) or {}
    pip_deps = [d for d in (spec.get("dependencies") or [])
                if isinstance(d, dict) and "pip" in d]
    if pip_deps:
        raise click.ClickException(
            f"{yaml_path.name} has a `pip:` subsection, but the CPU lock is rendered with "
            "`conda-lock --kind explicit` (conda-only) — pip deps would be silently dropped "
            "from the lock and missing from every provisioned env. Move those packages to "
            "conda dependencies (all of MuAgene's are on conda-forge/bioconda).")
    if not shutil.which("conda-lock"):
        raise click.ClickException(
            "conda-lock not found. Install dev deps:  pip install 'Processing-MuAgent[dev]'  "
            "(or: pip install conda-lock).")
    # Stamp the lock with the YAML's content hash; the env preflight compares this (not
    # mtimes — git doesn't preserve those) to detect a lock that drifted from the YAML.
    src_hash = hashlib.sha256(yaml_path.read_bytes()).hexdigest()
    for plat in platforms:
        click.echo(f"conda-lock --kind explicit -p {plat} -f {yaml_path}")
        try:
            subprocess.run(["conda-lock", "--kind", "explicit", "-f", str(yaml_path), "-p", plat],
                           cwd=str(work), check=True)
        except (subprocess.SubprocessError, OSError) as exc:
            raise click.ClickException(f"conda-lock failed for {plat}: {exc}")
        produced = work / f"conda-{plat}.lock"        # conda-lock's default explicit name
        dest = work / f"processing.{plat}.lock"        # the manifest's convention
        if produced.exists():
            produced.replace(dest)
        dest.write_text(f"# source-sha256: {src_hash}\n" + dest.read_text())
        click.echo(f"wrote {dest}")
    click.echo("Lockfile(s) regenerated. Commit them alongside processing.yaml.")


@main.command()
@click.argument("stage")
@click.argument("param_kv")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--rationale", default="User revision")
@click.option("--dry-run", is_flag=True,
              help="Preview the parameter change and exactly which artifacts would be "
                   "deleted (plus the current QC thresholds) — mutate nothing.")
def revise(stage: str, param_kv: str, config_path: str, rationale: str, dry_run: bool) -> None:
    """Update one parameter and reset the stage to awaiting_approval.

    PARAM_KV is key=value, e.g. s1_rna_qc.pct_counts_mt_max=10.0
    """
    stage = _canonical_stage(stage)
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    if "=" not in param_kv:
        raise click.ClickException("param_kv must be key=value")
    key, value = param_kv.split("=", 1)
    # Accept both the short form (min_counts_floor=500) and the full form
    # (s1_rna_qc.min_counts_floor=500). Without this normalisation the key is
    # stored bare in parameters.yaml and effective_params() — which looks for
    # "<stage>.<param>" — never finds it, so the revise has no effect on the
    # plan-review preview or on the real stage at runtime.
    if not key.startswith(f"{stage}."):
        key = f"{stage}.{key}"
    try:
        value_parsed = yaml.safe_load(value)
    except Exception:
        value_parsed = value
    if dry_run:
        _revise_dry_run(run_dir, paths, stage, key, value_parsed)
        return
    provenance.set_param(
        str(paths.parameters_yaml),
        key, value_parsed,
        source="user", confidence="high", rationale=rationale,
    )
    approval.mark_awaiting(run_dir, stage)
    log_event(run_dir, {"stage": stage, "event": "revised", "param": key, "value": value_parsed})
    click.echo(f"Revised {key} = {value_parsed!r}; {_display_stage(stage)} is awaiting_approval.")

    # A revise behaves differently depending on which checkpoint is active.
    if not approval.is_approved(run_dir, "plan_review"):
        # Plan-review checkpoint: the plan is not locked yet. The override only
        # tunes a proposed value, so re-render the plan deliverables (overlay) so
        # what the user reviews equals what will run. No QC stage has executed,
        # so there is nothing downstream to invalidate.
        regenerated = _regenerate_plan_deliverables(run_dir)
        if regenerated:
            log_event(run_dir, {"stage": "plan_review", "event": "plan_deliverables_regenerated",
                                "regenerated": regenerated})
            click.echo(
                "Regenerated plan deliverables so the review reflects this revise "
                f"(overlay): {', '.join(regenerated)}. plan_review stays awaiting_approval."
            )
        return

    # Post-approval (QC-review checkpoint): re-running a QC stage requires
    # clearing its stale downstream artifacts AND the post_qc_review gate
    # outputs, or Snakemake reports "Nothing to be done" and silently skips the
    # re-run. Do this deterministically here so the agent never hand-deletes.
    invalidated = _invalidate_qc_downstream(run_dir, stage)
    if invalidated:
        log_event(run_dir, {"stage": stage, "event": "qc_downstream_invalidated",
                            "deleted": invalidated})
        click.echo(
            f"Invalidated {len(invalidated)} stale downstream/gate artifact(s) so the "
            f"re-run regenerates them (incl. the post_qc_review gate). "
            f"Approve {_display_stage(stage)} (and s3_doublets if S1/S2 changed), then submit."
        )


def _stage_states(paths: RunPaths) -> list[tuple[str, str]]:
    return _sp.stage_states(paths)


@main.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--watch", is_flag=True,
              help="Poll until a review gate needs approval or manifest completes.")
@click.option("--interval", type=float, default=15.0,
              help="Poll interval in seconds when --watch is set.")
def status(config_path: str, watch: bool, interval: float) -> None:
    """Print per-step pipeline state (S1a–S8 + review gates). With --watch, polls until something changes."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)

    def _print(states: list[tuple[str, str, str]]) -> None:
        for label, task, st in states:
            click.echo(f"  {label:18s}  {task:30s}  {st}")

    if not watch:
        _print(_stage_states(paths))
        return

    sys.stdout.reconfigure(line_buffering=True)
    last: list[tuple[str, str, str]] | None = None
    while True:
        states = _stage_states(paths)
        if states != last:
            click.echo(f"--- {time.strftime('%Y-%m-%d %H:%M:%S')} ---")
            _print(states)
            last = states
        else:
            active = next((task for _, task, st in states if st == "in_progress"), "idle")
            click.echo(f"[{time.strftime('%H:%M:%S')}] {active}")
        if any(st == "failed" for _, _, st in states):
            click.echo("\n→ a step failed; inspect logs under "
                       f"{paths.snakemake_workdir}/.snakemake/slurm_logs/ "
                       "then fix and `submit` again (resume target is inferred).")
            return
        if any(st == "cancelled" for _, _, st in states):
            click.echo("\n→ a step was cancelled by the HPC monitor; see "
                       f"{paths.run_dir / 'internal' / 'hpc_monitor' / 'latest_report.md'} "
                       "for the confirmed-dead reason and investigation evidence. "
                       "Re-`submit` to resume.")
            return
        if any(st == "awaiting_approval" for _, _, st in states):
            click.echo("\n→ a review gate is awaiting approval; review deliverables and run "
                       "`Processing-MuAgent approve <stage>` (e.g. qc_review).")
            return
        if paths.run_manifest_json.exists():
            click.echo("\n→ run_manifest.json present; pipeline complete.")
            return
        time.sleep(max(2.0, interval))


@main.command(name="hpc-status")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def hpc_status(config_path: str) -> None:
    """Report HPC job health, monitor findings, and per-step pipeline state (one-shot).

    This is Processing-MuAgent's single window onto the Execution-MuAgent supervision
    daemon, which is the sole monitor. It reads only structured JSON
    (latest_snapshot.json + latest_submission.json) — health, silence/tolerance,
    findings, kill_action, and supervisor liveness — and prints once, then exits.

    There is no poll loop here: the daemon does the monitoring. This command drives
    the report-and-repoll rule — after `submit`, report this status, then (while the
    job is still running) re-poll on a non-blocking scheduled wakeup after the seconds
    printed on the `Next check:` line, until monitor.pid is removed or a review gate is
    awaiting approval (`Gate signal present`). Report to the user only when the `State:`
    fingerprint changes.
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    submission = _sp.load_latest_hpc_submission(paths)
    snapshot = _sp.load_hpc_monitor_state(paths) or {}

    def _supervisor_status() -> str:
        pid_path = paths.run_dir / "internal" / "hpc_monitor" / "monitor.pid"
        if not pid_path.exists():
            return "not started"
        try:
            pid = int(pid_path.read_text().strip())
            os.kill(pid, 0)
            return f"alive (PID {pid})"
        except (ValueError, OSError):
            return "not running"

    def _print_hpc_header() -> None:
        if submission is None:
            click.echo("  (no HPC submission registered)")
            return
        job_id = submission.get("job_id", "?")
        stage = submission.get("target", "?").removesuffix("_execute")
        submitted_at = submission.get("submitted_at", "?")
        sched = (snapshot.get("scheduler") or {}) if snapshot else {}
        sched_state = sched.get("state") or "unknown"
        elapsed = sched.get("elapsed") or "?"
        timelimit = sched.get("timelimit") or "?"
        ms = (snapshot.get("monitor_state") or {}) if snapshot else {}
        health = ms.get("health", "unknown")
        silence = ms.get("silence_intervals", "?")
        tolerance = ms.get("tolerance_n", "?")
        sup_status = _supervisor_status()
        sup_offline = "not running" in sup_status
        displayed_health = "stale (supervisor offline)" if sup_offline else health
        click.echo(f"  Stage: {stage}  Job: {job_id}  Submitted: {submitted_at}")
        click.echo(f"  Scheduler: {sched_state}  Elapsed: {elapsed} / {timelimit}")
        click.echo(f"  Health: {displayed_health}  (silence {silence}/{tolerance} intervals)")
        click.echo(f"  Supervisor: {sup_status}")
        reason = ms.get("confirmed_dead_reason")
        if reason:
            click.echo(f"  Confirmed-dead reason: {reason}")
        if "not running" in sup_status:
            active_states = {"running", "pending", "r", "q", "cf"}
            if any(s in sched_state.lower() for s in active_states):
                click.echo(
                    "\n  WARNING: Supervision offline — stalled/failed jobs will NOT be "
                    "auto-cancelled. Restart with:\n"
                    f"    Processing-MuAgent supervisor-restart --config {paths.run_yaml}"
                )

    def _print_findings() -> None:
        findings = _sp.load_hpc_findings(paths)
        if not findings:
            return
        click.echo("")
        click.echo("--- Monitor findings (latest check) ---")
        for f in findings:
            severity = str(f.get("severity", "info")).upper()
            click.echo(f"  [{severity}] {f.get('code', '?')}: {f.get('message', '')}")

    def _print(states: list[tuple[str, str, str]]) -> None:
        click.echo("")
        click.echo("--- HPC monitor ---")
        _print_hpc_header()
        _print_findings()
        click.echo("")
        click.echo("--- Pipeline state ---")
        for label, task, st in states:
            click.echo(f"  {label:18s}  {task:30s}  {st}")

    states = _stage_states(paths)
    _print(states)

    # --- report-and-repoll: deterministic fingerprint + the single next-action ---
    state_set = {st for _, _, st in states}
    gate = "awaiting_approval" in state_set
    cancelled = "cancelled" in state_set
    failed = "failed" in state_set
    complete = paths.run_manifest_json.exists()
    pid_present = (paths.run_dir / "internal" / "hpc_monitor" / "monitor.pid").exists()
    sup_alive = pid_present and "alive" in _supervisor_status()

    sched_state = ((snapshot.get("scheduler") or {}).get("state")) or "unknown"
    health = (snapshot.get("monitor_state") or {}).get("health", "unknown")
    n_findings = len(_sp.load_hpc_findings(paths) or [])
    # Stable single line the report-and-repoll rule diffs to decide whether to re-report.
    click.echo("")
    click.echo(f"State: {sched_state}/{health}/sup={'alive' if sup_alive else 'offline'}/"
               f"gate={'awaiting_approval' if gate else 'none'}/findings={n_findings}")

    # Informative guidance for terminal / gate states.
    if cancelled:
        click.echo("\n→ a step was cancelled by the HPC monitor (see the kill_action / "
                   "confirmed-dead reason above). Fix the cause and re-`submit` to resume.")
    elif failed:
        click.echo("\n→ a step failed; inspect logs under "
                   f"{paths.snakemake_workdir}/.snakemake/slurm_logs/ "
                   "then fix and `submit` again (resume target is inferred).")
    elif gate:
        click.echo("\n→ a review gate is awaiting approval; review deliverables and run "
                   "`Processing-MuAgent approve <stage>` (e.g. qc_review).")
    elif complete:
        click.echo("\n→ run_manifest.json present; pipeline complete.")

    # The single next-action token the report-and-repoll rule keys on.
    if gate or complete:
        click.echo("→ Gate signal present — drive the next checkpoint now (no re-poll).")
    elif cancelled or failed:
        click.echo("→ Run halted — no re-poll; resolve the above, then `submit` again.")
    elif sup_alive:
        repoll_s = int(snapshot.get("next_recheck_after_s") or hpc.DEFAULT_REPOLL_AFTER_S)
        interval_s = int(snapshot.get("interval_s") or 270)
        click.echo(f"Next check: re-poll via scheduled wakeup in ~{repoll_s}s "
                   f"(daemon interval {interval_s}s + {repoll_s - interval_s}s buffer)")
    else:
        click.echo("→ Supervisor offline with no gate armed — re-check pipeline state / logs "
                   "before re-polling, or `supervisor-restart` if the cluster job is still active.")


@main.command(name="supervisor-restart")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--kill-existing/--no-kill-existing", default=True, show_default=True,
              help="Kill any running supervisor before starting the new one.")
def supervisor_restart(config_path: str, kill_existing: bool) -> None:
    """Restart the background supervisor daemon without resubmitting the cluster job.

    Use when the supervisor process died mid-run (crash, OOM, site reboot) but the
    cluster job is still active. Reads latest_submission.json and re-invokes
    resume-monitor as a new daemon (no resubmit).

    The supervisor is the kill-on-hang safety layer. Restarting it restores stall
    detection and auto-cancel protection for the running job.
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    sub_path = paths.run_dir / "internal" / "hpc_monitor" / "latest_submission.json"
    if not sub_path.exists():
        raise click.ClickException(
            "No submission recorded for this run. Use `submit` to start a job."
        )
    if kill_existing:
        killed = hpc.kill_existing_supervisor(run_dir)
        if killed:
            click.echo("Stopped existing supervisor.")
    env = hpc._execution_muagent_env()
    if env is None:
        raise click.ClickException(
            "Execution-MuAgent not found. Install it: pip install -e Execution-MuAgent/"
        )
    cmd = [
        sys.executable, "-m", "execution_muagent.cli", "resume-monitor",
        "--run-dir", str(run_dir),
    ]
    result = hpc.start_supervisor_daemon(run_dir, cmd, env)
    if result is None:
        raise click.ClickException("Failed to start supervisor daemon.")
    pid = result["pid"]
    log = result["log"]
    log_event(run_dir, {
        "stage": "supervisor_restart", "event": "restarted",
        "supervisor_pid": pid, "supervisor_log": log,
    })
    click.echo(f"Supervisor restarted (PID {pid}), logging to {log}")
    click.echo(f"Report status: Processing-MuAgent hpc-status --config {paths.run_yaml}")


@main.command(name="plan-review")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--intro", "intro_text", default=None,
              help="Introductory paragraph to prepend before the Summary section.")
@click.option("--intro-context", "intro_context_only", is_flag=True, default=False,
              help="Print the intro context JSON and exit without writing plan_review.md.")
def plan_review_cmd(config_path: str, intro_text: str | None, intro_context_only: bool) -> None:
    """Render and write the merged plan-review markdown (summary + appendix).

    Also writes per-stage job spec YAMLs to internal/specs/ so Execution-MuAgent
    can read science intent, resource hints, and progress_timeout_hint per stage.

    This command is a *renderer*: it requires the planning compute (P1 → S0 → P2)
    to have finished and produced preprocessing_plan.json. Calling it before that
    would emit placeholder deliverables and a false awaiting_approval signal.
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    if intro_context_only:
        import json as _json
        click.echo(_json.dumps(_pr.build_intro_context(run_dir), indent=2))
        return
    # Guard against rendering before planning compute has produced the plan.
    missing = []
    if not paths.preprocessing_plan.exists():
        missing.append(str(paths.preprocessing_plan))
    if not paths.validation_report.exists():
        missing.append(str(paths.validation_report))
    if missing:
        raise click.ClickException(
            "Cannot render plan review — S0 ingest has not finished yet.\n"
            f"Missing: {', '.join(missing)}\n"
            "Wait for the planning job (target plan_review_propose) to complete, "
            "then re-run this command."
        )
    text = _pr.render_merged_markdown(run_dir, intro=intro_text)
    click.echo(text)
    out = _pr.write_summary(run_dir, intro=intro_text)
    click.echo(f"\nWritten: {out}")
    html_out = _pr.write_plan_summary_html(run_dir, intro=intro_text)
    click.echo(f"Written: {html_out}")
    # Write per-stage specs; read workflow_branch from plan if available.
    try:
        import json
        plan_path = RunPaths(run_dir).preprocessing_plan
        branch = "paired"
        if plan_path.exists():
            branch = json.loads(plan_path.read_text()).get("workflow_branch", "paired")
        written = _specs.write_stage_specs(run_dir, branch)
        if written:
            click.echo(f"Wrote {len(written)} stage metadata file(s) to {RunPaths(run_dir).stage_meta_dir}/")
    except Exception:
        pass  # spec writing is best-effort; never block plan-review
    # Arm the plan_review gate only when the plan actually exists. The primary
    # gate-arming path is the plan_review_propose Snakemake rule; this CLI path
    # is a re-render convenience after planning compute has finished.
    approval.mark_awaiting(run_dir, "plan_review")


def _unlock_snakemake(run_dir: Path, config_path: Path) -> None:
    paths = RunPaths(run_dir)
    env = os.environ.copy()
    env.setdefault("PYTHONPATH", str(PACKAGE_DIR))
    env.setdefault("PMA_REPO_ROOT", str(PACKAGE_DIR))
    paths.snakemake_workdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "-m", "snakemake",
        "-s", str(SNAKEFILE),
        "--directory", str(paths.snakemake_workdir),
        "--unlock",
        "--configfile", str(config_path),
    ]
    click.echo(f"$ {' '.join(cmd)}")
    result = subprocess.run(cmd, env=env, cwd=str(PACKAGE_DIR))
    if result.returncode != 0:
        raise click.ClickException(f"snakemake --unlock exited with {result.returncode}")


@main.command(name="marker-gene-check")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option(
    "--force-tsne",
    is_flag=True,
    default=False,
    help="Recompute t-SNE even when a valid cache exists.",
)
@click.option(
    "--plot-only",
    is_flag=True,
    default=False,
    help="Write the figure only; do not refresh QC review reports.",
)
@click.argument("genes", nargs=-1, required=True)
def marker_gene_check_cmd(
    config_path: str,
    force_tsne: bool,
    plot_only: bool,
    genes: tuple[str, ...],
) -> None:
    """Generate before/after marker gene expression plots.

    GENES is one or more gene symbols, e.g. CD3E CD20 EPCAM (matched case-insensitively).

    Uses a cached t-SNE embedding when the cell set is unchanged. By default, QC review
    reports are refreshed automatically after plotting. Pass ``--plot-only`` to skip that.
    """
    from .stages import s1a_ambient as _s1a
    run_dir = _resolve_run_dir(config_path)

    if not genes:
        raise click.UsageError("Provide at least one gene symbol.")

    gene_list = list(genes)
    click.echo(f"Checking marker genes: {', '.join(gene_list)}")
    result = _s1a.run_marker_gene_check(
        run_dir, gene_list, force_tsne=force_tsne, refresh_qc=not plot_only,
    )
    if result["found"]:
        click.echo(f"Plotted: {', '.join(result['found'])}")
    else:
        click.echo("No marker genes found in matrix; figure not written.")
    if result["missing"]:
        click.echo(f"Not found in data: {', '.join(result['missing'])}")
    if result["found"]:
        if plot_only:
            click.echo("Figure written (--plot-only: QC reports unchanged).")
        else:
            click.echo("QC reports refreshed.")


@main.command(name="unlock")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def unlock_cmd(config_path: str) -> None:
    """Remove stale Snakemake locks for a run after confirming no active process."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    locks = hpc.snakemake_lock_files(paths.snakemake_workdir)
    if not locks:
        click.echo(f"No Snakemake locks found under {paths.snakemake_workdir}.")
        return
    active = hpc.snakemake_processes_for_workdir(paths.snakemake_workdir)
    if active:
        detail = "\n".join(f"  pid {pid}: {args}" for pid, args in active)
        raise click.ClickException(
            "Refusing to unlock while a local Snakemake process references this workdir:\n"
            f"{detail}"
        )
    _unlock_snakemake(run_dir, Path(config_path))
    click.echo(f"Unlocked {paths.snakemake_workdir}")


def _enforce_context_gate(paths: RunPaths, no_context: bool) -> None:
    """Phase 1 biological-context gate (MANDATORY before any preprocessing).

    Shared by `run` (local execution) and `submit` (cluster execution) so the gate
    is enforced identically regardless of which entry point starts the pipeline.
    Raises if the report is still the blank template and the caller did not pass
    --no-context.
    """
    report_path = paths.biological_context_md
    report_path.parent.mkdir(parents=True, exist_ok=True)
    if not report_path.exists():
        _ctx.write_template(report_path)
    if _ctx.is_unfilled_template(report_path) and not no_context:
        raise click.ClickException(
            "Biological Context Report at "
            f"{report_path}\nis empty (template only). Per Phase 1 policy, preprocessing "
            "cannot proceed until the user provides biological context.\n\n"
            "Choose one of:\n"
            "  1. Paste context into the report file (fields: Organism, Tissue / sample, "
            "Assay, DOI(s) optional, Notes optional) and re-run.\n"
            "  2. Supply a report document (.docx/.pdf/.md/.txt) path in run.yaml under "
            "'biological_context_path' and re-run.\n"
            "  3. Explicitly proceed without context by adding --no-context to this command; "
            "the subagent will mark user-declared fields as status=missing and rely on "
            "file inputs + inference.\n"
        )
    if _ctx.is_unfilled_template(report_path) and no_context:
        click.echo("Proceeding WITHOUT biological context (--no-context set). User-declared "
                   "fields will be marked status=missing.", err=True)


def _enforce_execution_mode_gate(run_dir: Path, paths: RunPaths) -> None:
    """Require an explicit, user-confirmed execution mode before launching compute.

    System requirement: Processing-MuAgent must ALWAYS confirm execution mode
    (local vs HPC) with the user before running ANY compute job — not only at S0.
    This is enforced unconditionally on every `run` and `submit`, so it also
    covers resume submissions (S1+) and runs whose config never recorded an
    execution mode. Mirrors `_enforce_context_gate`.

    Idempotent: once `execution.user_confirmed` is true the gate passes instantly
    on every subsequent call, so re-checking has no cost.
    """
    params = str(paths.parameters_yaml)
    mode = provenance.get_value(params, "execution.mode", None)
    confirmed = provenance.get_value(params, "execution.user_confirmed", False)
    if mode is None:
        raise click.ClickException(
            "Execution mode is not set for this run. Per system policy, you MUST "
            "confirm local vs HPC with the user before any compute job runs "
            "(this applies to resume sessions too, not only the first S0 ingest).\n"
            "  1. Ask the user: run locally on this machine, or submit to an HPC "
            "cluster (SLURM)?\n"
            "  2. For HPC, probe the login node first: Processing-MuAgent hpc-info\n"
            "  3. Record the user's explicit choice:\n"
            f"     Processing-MuAgent configure-execution --config {paths.run_yaml} "
            "--mode local --confirmed-by-user\n"
            "     (or --mode slurm with partition + account, "
            "plus --confirmed-by-user)"
        )
    if not confirmed:
        raise click.ClickException(
            f"Execution mode is set to {mode!r} but was not confirmed by the user. "
            "Per system policy, confirm local vs HPC with the user before launching "
            "any compute job, then re-run:\n"
            f"  Processing-MuAgent configure-execution --config {paths.run_yaml} "
            f"--mode {mode} --confirmed-by-user"
        )


@main.command(name="run")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--auto-approve", is_flag=True, help="Auto-approve every checkpoint (noninteractive).")
@click.option("--auto-approve-except", "auto_except", multiple=True,
              help="With --auto-approve, do NOT pre-seed the given stage(s). Repeatable. "
                   "Example: --auto-approve-except qc_review")
@click.option("--no-context", is_flag=True, help="Explicit user choice to proceed without biological context; fields marked status=missing.")
@click.option("--marker-genes", "marker_genes_ack",
              type=click.Choice(["defer", "skip"]), default=None,
              help="With --auto-approve: record an explicit marker-gene decision so "
                   "plan_review can be seeded. 'defer' = check at QC review; "
                   "'skip' = decline. Provide actual genes via `revise` instead to run the check.")
@click.option("--target", default="all")
def run_pipeline(config_path: str, auto_approve: bool, auto_except: tuple[str, ...],
                 no_context: bool, marker_genes_ack: str | None, target: str) -> None:
    """Run the DAG LOCALLY. With --auto-approve, checkpoints are unblocked automatically.

    `run` is local-only: it executes on this machine (local mode) or runs the
    login-node localrules (propose / planning / manifest). All cluster job
    submission and monitoring is owned by Execution-MuAgent via `submit` — there
    is no `run --executor slurm` path.

    Use --auto-approve-except <stage> to keep specific gates honoured (e.g.
    qc_review in headless HPC mode).
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)

    _enforce_execution_mode_gate(run_dir, paths)
    mode = provenance.get_value(str(paths.parameters_yaml), "execution.mode", None)
    if mode == "slurm":
        raise click.ClickException(
            f"execution.mode is {mode!r}, but `run` is local-only. Heavy stages "
            "(starting with S0 ingest) must run on a compute node, never the login "
            "node. Submit instead:\n"
            f"  source {paths.hpc_env_sh}\n"
            f"  Processing-MuAgent submit --config {paths.run_yaml} "
            f"--executor {mode} --target {target}"
        )

    _enforce_context_gate(paths, no_context)

    auto_except = tuple(_canonical_stage(s) for s in auto_except)
    if auto_approve:
        # Pre-seed approval sentinels so snakemake can run the DAG end-to-end in a
        # single invocation; --auto-approve-except keeps the listed stages gated.
        kept = set(auto_except)
        _apply_marker_gene_ack(run_dir, marker_genes_ack)
        _seed_approvals(run_dir, HUMAN_CHECKPOINT_STAGES, note="auto-approved", kept=kept)
        if kept:
            click.echo(f"Auto-approved all stages except: "
                       f"{sorted(_display_stage(s) for s in kept)}. "
                       "Snakemake will stop at those gates.")
    _snakemake(["--configfile", str(paths.run_yaml), target], run_dir)


def _infer_submit_target(run_dir: Path) -> str:
    """Pick the Snakemake target from the first incomplete pipeline step."""
    return _sp.infer_resume_target(run_dir)


def _missing_approvals(run_dir: Path, stages: tuple[str, ...]) -> list[str]:
    paths = RunPaths(run_dir)
    return [stage for stage in stages if not paths.approved_sentinel(stage).exists()]


def _seed_approvals(
    run_dir: Path,
    stages: tuple[str, ...],
    *,
    note: str,
    kept: set[str] | None = None,
) -> list[str]:
    kept = kept or set()
    seeded: list[str] = []
    protected = False  # any gate left in the approved state (freshly seeded OR already)
    for stage in stages:
        if stage in kept:
            continue
        if approval.is_approved(run_dir, stage):
            # Already approved on a prior call — do NOT re-stamp. approval.approve
            # rewrites the sentinel and bumps its mtime; since every QC/downstream
            # execute rule declares <gate>.approved as an input, a fresh mtime is an
            # "Updated input files" trigger that needlessly re-runs the whole approved
            # upstream chain on the next submit. The gate is already honoured; leave
            # the original sentinel (and its mtime) untouched.
            protected = True
            continue
        if stage == "plan_review" and _pr.marker_gene_decision_pending(run_dir):
            raise click.ClickException(
                "Cannot auto-approve plan_review: " + _MARKER_GENE_GATE_MSG
                + "\nFor an unattended batch, pass --marker-genes defer|skip.")
        approval.approve(run_dir, stage, note=note)
        seeded.append(stage)
        protected = True
        if stage == "post_qc_review":
            deleted = _cleanup_qc_intermediates(run_dir)
            if deleted:
                log_event(run_dir, {"stage": "post_qc_review", "event": "qc_cleanup",
                                     "deleted": deleted})
    # Set the revoke-protection flag whenever a gate is in the approved state, not
    # only when one was freshly seeded — otherwise a submit that re-enters an
    # already-approved phase would let the propose rules revoke those approvals.
    if protected:
        os.environ["PMA_AUTO_APPROVE"] = "1"
    return seeded


def _prepare_submit_approvals(
    run_dir: Path,
    target: str,
    *,
    inferred_target: bool,
    auto_approve: bool,
    auto_except: tuple[str, ...],
) -> list[str]:
    """Seed internal phase approvals or fail fast for explicit unsafe targets."""
    internal: tuple[str, ...] = ()
    human = _sp.required_human_approvals(target)
    missing_human = _missing_approvals(run_dir, human)
    if missing_human:
        raise click.ClickException(
            f"Target {target!r} requires human approval sentinel(s): "
            f"{', '.join(_display_stage(s) for s in missing_human)}. "
            "Review/approve these gates before submitting."
        )

    if auto_approve:
        return []

    if inferred_target:
        kept = set(auto_except)
        seeded = _seed_approvals(
            run_dir,
            tuple(s for s in internal if s not in _sp.required_human_approvals("all")),
            note=f"phase-auto-approved for {target}",
            kept=kept,
        )
        return seeded

    missing_internal = _missing_approvals(run_dir, internal)
    if missing_internal:
        raise click.ClickException(
            f"Explicit target {target!r} requires internal approval sentinel(s): "
            f"{', '.join(missing_internal)}. Use --auto-approve for an unattended "
            "batch, approve those stages, or omit --target so submit can infer and "
            "prepare the current phase."
        )
    return []


@main.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--executor", type=CLUSTER_EXECUTOR_CHOICE, required=True,
              help="Scheduler to submit the head-job to (slurm). "
                   "For local foreground runs use `run` (which is local-only).")
@click.option("--target", default=None,
              help="Override the Snakemake target. Omit to auto-infer the first "
                   "incomplete step (e.g. plan_review_propose for planning, "
                   "post_qc_review_propose, all).")
@click.option("--no-context", is_flag=True,
              help="Explicit user choice to proceed without biological context (planning "
                   "submissions only); fields marked status=missing.")
@click.option("--auto-approve", is_flag=True,
              help="Pre-seed all checkpoint sentinels; head-job runs unattended end-to-end.")
@click.option("--auto-approve-except", "auto_except", multiple=True,
              help="With --auto-approve, keep these gates honoured. Repeatable.")
@click.option("--marker-genes", "marker_genes_ack",
              type=click.Choice(["defer", "skip"]), default=None,
              help="With --auto-approve: record an explicit marker-gene decision so "
                   "plan_review can be seeded. 'defer' = check at QC review; "
                   "'skip' = decline. Provide actual genes via `revise` instead to run the check.")
@click.option("--output", "output_log", type=click.Path(), default=None,
              help="Scheduler output-log path for the head-job (optional).")
@click.option("--unlock-stale-locks", is_flag=True,
              help="If Snakemake locks exist and no local process owns this workdir, "
                   "run snakemake --unlock before submitting.")
@click.option("--watch/--no-watch", default=True,
              help="Start a background supervisor daemon that monitors the cluster job and "
                   "cancels it if it hangs (default: on). The daemon survives SSH disconnect "
                   "unless the site uses KillUserProcesses=yes — use tmux/screen there. "
                   "Returns after job submission is confirmed (≤90 s). "
                   "--no-watch: submit only, NO supervisor daemon started — no stall "
                   "detection, no auto-cancel.")
def submit(config_path: str, executor: str, target: str | None, no_context: bool,
           auto_approve: bool, auto_except: tuple[str, ...],
           marker_genes_ack: str | None,
           output_log: str | None, unlock_stale_locks: bool,
           watch: bool) -> None:
    """Submit the snakemake runner as a SLURM head-job.

    This is the ONLY cluster-execution path: Processing-MuAgent prepares the
    head-job spec + site.config and Execution-MuAgent owns submission and
    monitoring (kill-on-hang, hpc-status). The planning phase targets
    ``plan_review_propose`` (auto-inferred), which pulls P1 → S0 → P2 as
    Snakemake dependencies and arms the gate at the end of a single head-job.

    Execution-MuAgent is a hard dependency for cluster submission — it renders the
    submission script, submits the head-job, and owns monitoring. If Execution-MuAgent
    is unavailable, this command fails loudly: there is no manual-submission path.

    The head-job runs on a compute node, activates the project conda env, and
    invokes snakemake with the cluster profile. Snakemake then submits per-stage
    child jobs. The head-job exits when the DAG completes or stops at a missing
    approval gate.

    Typical headless workflow on HPC:

        # Run planning interactively (Phase A), then submit the heavy middle:
        Processing-MuAgent submit --config $CFG --executor slurm \\
                --auto-approve --auto-approve-except post_qc_review

        # After QC review, approve and resume (target auto-inferred); the run
        # proceeds through clustering and UMAP to the final outputs with no
        # further checkpoint:
        Processing-MuAgent approve post_qc_review --config $CFG
        Processing-MuAgent submit --config $CFG --executor slurm
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)

    # System requirement: confirm execution mode with the user before launching ANY
    # cluster job. Fires before approval seeding, target inference, and the
    # site.config check — so resume submissions (S1+) are gated too, not only S0.
    _enforce_execution_mode_gate(run_dir, paths)

    auto_except = tuple(_canonical_stage(s) for s in auto_except)
    if auto_approve:
        kept = set(auto_except)
        _apply_marker_gene_ack(run_dir, marker_genes_ack)
        _seed_approvals(run_dir, HUMAN_CHECKPOINT_STAGES, note="auto-approved (submit)", kept=kept)
        if kept:
            click.echo(f"Auto-approved all stages except: "
                       f"{sorted(_display_stage(s) for s in kept)}.")
        # Tell the head-job's propose rules not to revoke pre-seeded approvals.
        os.environ["PMA_AUTO_APPROVE"] = "1"

    inferred_target = target is None
    resolved_target = target if target is not None else _infer_submit_target(run_dir)

    # Phase 1 biological-context gate — enforced for planning-phase submissions
    # exactly as `run` does. Resume submissions (S1+) skip it: context was
    # already validated when planning ran. Covers both the canonical auto-inferred
    # target (plan_review_propose) and legacy explicit --target s0_ingest_execute.
    if resolved_target in {"s0_ingest_execute", "plan_review_propose"}:
        _enforce_context_gate(paths, no_context)

    phase_seeded = _prepare_submit_approvals(
        run_dir,
        resolved_target,
        inferred_target=inferred_target,
        auto_approve=auto_approve,
        auto_except=auto_except,
    )

    locks = hpc.snakemake_lock_files(paths.snakemake_workdir)
    if locks:
        active = hpc.snakemake_processes_for_workdir(paths.snakemake_workdir)
        if active:
            detail = "\n".join(f"  pid {pid}: {args}" for pid, args in active)
            raise click.ClickException(
                "Snakemake locks exist and a local Snakemake process still references "
                f"{paths.snakemake_workdir}:\n{detail}"
            )
        lock_list = ", ".join(str(p) for p in locks)
        if not unlock_stale_locks:
            raise click.ClickException(
                "Snakemake lock files already exist for this run, so submitting now "
                "would fail with LockException.\n"
                f"Locks: {lock_list}\n"
                "If no scheduler head/child jobs for this run are active, recover with:\n"
                f"  Processing-MuAgent unlock --config {paths.run_yaml}\n"
                "or resubmit with `--unlock-stale-locks`."
            )
        click.echo(f"Unlocking stale Snakemake locks: {lock_list}")
        _unlock_snakemake(run_dir, paths.run_yaml)

    # Archive the prior run's Snakemake logs so a resubmit's PENDING window is not
    # misread as the previous run's failure by `hpc-status` (stage state is derived
    # from the newest per-rule + main logs). No-op on a first submit.
    archived = hpc.archive_prior_run_logs(paths.snakemake_workdir)
    if archived is not None:
        click.echo(f"Archived previous run logs → {archived}")

    out_path = Path(output_log) if output_log else hpc.head_job_log_path(executor)

    if not paths.site_config.exists():
        raise click.ClickException(
            f"site.config not found at {paths.site_config}. "
            "Run `Processing-MuAgent configure-execution --mode slurm ...` first."
        )

    # Write the head-job spec so Execution-MuAgent can render + submit it.
    head_spec_path = _specs.write_head_job_spec(run_dir, resolved_target)

    if watch:
        click.echo("Starting supervision daemon (background)...")

    ea_result = hpc.submit_via_execution_muagent(
        head_spec_path,
        paths.site_config,
        run_dir,
        resolved_target,
        watch=watch,
        kill_on_hang=True,
    )
    if ea_result is None:
        raise click.ClickException(
            "Execution-MuAgent is required for cluster submission but is not available "
            "or returned an error.\n"
            "  Install it:  pip install -e Execution-MuAgent/\n"
            "  Then re-run `Processing-MuAgent submit`."
        )

    import re as _re
    if watch:
        pid = ea_result.get("pid", "?")
        mon_log = ea_result.get("log", "")
        # ea_result["job_id"] is the entry the supervisor appended for THIS
        # submission (or None on timeout). Do NOT fall back to the last manifest
        # entry — that can be a stale job_id from a previous head-job.
        job_id = ea_result.get("job_id")
        if not job_id:
            click.echo(
                "Warning: job ID not yet in execution_manifest.jsonl "
                "(scheduler slow or NFS lag). Check `hpc-status` to confirm submission.",
                err=True,
            )
            job_id = "unknown"
        submitted_log_path = hpc.submitted_log_path(executor, out_path, job_id)
        log_event(run_dir, {
            "stage": "submit", "event": "head_job_submitted",
            "executor": executor, "target": resolved_target,
            "job_id": job_id, "auto_approve": auto_approve,
            "kept_gates": sorted(set(auto_except)),
            "phase_auto_approved": phase_seeded,
            "via_execution_agent": True,
            "supervisor_pid": pid, "supervisor_log": mon_log,
            "head_job_log": str(submitted_log_path),
        })
        click.echo(f"Submitted {executor} head-job: {job_id}")
        click.echo(f"  config:     {paths.run_yaml}")
        click.echo(f"  target:     {resolved_target}")
        if phase_seeded:
            click.echo(f"  phase-auto-approved: {', '.join(phase_seeded)}")
        click.echo(f"  log:        {submitted_log_path}")
        click.echo(f"  supervisor: PID {pid}, log → {mon_log}")
        click.echo(
            "\nThe supervisor daemon is the sole monitor (kill-on-hang) and runs in "
            "the background; it survives SSH disconnect (unless the site uses "
            "KillUserProcesses=yes). Do not run a watch loop — report its status on "
            "demand and act when it signals (job terminal, or a review gate awaiting):\n"
            f"  Report status: Processing-MuAgent hpc-status --config {paths.run_yaml}\n"
            "The daemon signals completion by removing internal/hpc_monitor/monitor.pid; "
            "a review gate shows as 'awaiting_approval'."
        )
    else:
        m = _re.search(r"(?:head-job|job[_-]id)[:\s]+(\S+)", ea_result.get("stdout", ""))
        job_id = m.group(1).strip() if m else ea_result.get("stdout", "").strip().splitlines()[-1]
        submitted_log_path = hpc.submitted_log_path(executor, out_path, job_id)
        log_event(run_dir, {
            "stage": "submit", "event": "head_job_submitted",
            "executor": executor, "target": resolved_target,
            "job_id": job_id, "auto_approve": auto_approve,
            "kept_gates": sorted(set(auto_except)),
            "phase_auto_approved": phase_seeded,
            "via_execution_agent": True,
            "head_job_log": str(submitted_log_path),
        })
        click.echo(f"Submitted {executor} head-job: {job_id} (via Execution-MuAgent)")
        click.echo(f"  config:  {paths.run_yaml}")
        click.echo(f"  target:  {resolved_target}")
        if phase_seeded:
            click.echo(f"  phase-auto-approved: {', '.join(phase_seeded)}")
        click.echo(f"  log:     {submitted_log_path}")
        click.echo(
            "\nNote: --no-watch — no supervisor daemon started. "
            "Stalled/hung jobs will NOT be auto-cancelled.\n"
            f"Poll progress: Processing-MuAgent status --watch --config {paths.run_yaml}"
        )


def _snakemake(args: list[str], run_dir: Path) -> None:
    """Invoke snakemake LOCALLY with --cores 1 for reproducibility.

    `run` and `propose` are local-only entry points. All cluster execution is
    owned by Execution-MuAgent and reached via `submit` (which renders + submits
    a supervised head-job) — never through this helper. The head-job's own
    snakemake invocation (in launch_runner.sh) attaches the cluster profile; this
    helper does not.

    Expected args shape from callers: ["--configfile", <path>, <target>].
    """
    env = os.environ.copy()
    env.setdefault("PYTHONPATH", str(PACKAGE_DIR))
    env.setdefault("PMA_REPO_ROOT", str(PACKAGE_DIR))
    paths = RunPaths(run_dir)
    paths.snakemake_workdir.mkdir(parents=True, exist_ok=True)
    env.setdefault("XDG_CACHE_HOME", str(paths.snakemake_workdir / "cache"))
    # Single-thread for reproducibility (UMAP / numba) — unchanged on local;
    # cluster jobs inherit these unless the user overrides in their shell.
    env.setdefault("NUMBA_NUM_THREADS", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("PYTHONHASHSEED", "0")
    if os.environ.get("PMA_AUTO_APPROVE"):
        env["PMA_AUTO_APPROVE"] = os.environ["PMA_AUTO_APPROVE"]

    configfile_path = None
    targets: list[str] = []
    rest: list[str] = []
    it = iter(args)
    for a in it:
        if a == "--configfile":
            configfile_path = next(it, None)
        elif a.startswith("-"):
            rest.append(a)
        else:
            targets.append(a)

    cmd = [
        sys.executable, "-m", "snakemake",
        "-s", str(SNAKEFILE),
        "--directory", str(paths.snakemake_workdir),
        # Rerun only on mtime/missing-output, not params/code/input-set/software-env.
        # This pipeline forces reruns by explicit artifact deletion (executor revise
        # -> _invalidate_qc_downstream), so content/input-set triggers only cause
        # spurious reruns. Mirrors `rerun-triggers: [mtime]` in the cluster profiles.
        "--rerun-triggers", "mtime",
        "--rerun-incomplete", *targets, *rest,
        "--cores", "1",
    ]

    if configfile_path:
        cmd += ["--configfile", configfile_path]
    click.echo(f"$ {' '.join(cmd)}")
    r = subprocess.run(cmd, env=env, cwd=str(PACKAGE_DIR))
    if r.returncode != 0:
        raise click.ClickException(f"snakemake exited with {r.returncode}")


if __name__ == "__main__":
    main()
