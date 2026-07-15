"""Parameter revision policy — the QC-invalidation cascade + revise previews.

A `revise <stage> key=value` mutates parameters.yaml and, at the post_qc_review gate,
deterministically deletes the revised stage's outputs and everything strictly downstream
through S3 (plus the gate outputs) so a re-run actually re-executes. At plan_review nothing
has run, so a revise only re-renders the plan deliverables (no deletion). The dry-run path
previews exactly what a real revise would remove without mutating anything.
"""
from __future__ import annotations

from pathlib import Path

import click

from . import approval, provenance
from .cleanup import S1A_REGEN_CACHES
from .log import log_event
from .run_paths import RunPaths


# Per-stage QC artifacts that become stale when that stage is re-run, listed as
# (artifact_stage, filename). Revising a QC stage invalidates its own outputs and
# everything strictly downstream of it through S3 (the DAG is
# s1a_ambient -> s1_rna_qc -> s3_doublets and s2_atac_qc -> s3_doublets). The
# expensive chr-normalized fragment cache (atac_fragments_cbf_chrnorm.tsv.gz*) is
# intentionally NOT listed — it is reused across re-runs, and only deleted once QC
# is approved (by cleanup.cleanup_qc_intermediates, when no further re-run can occur).
_S1A_QC_ARTIFACTS = [
    ("s1a_ambient", "rna_decontaminated.h5ad"),  # untracked working file (read by S1 by path)
    ("s1a_ambient", "summary.json"),             # stage-done marker — delete to force re-run
] + S1A_REGEN_CACHES
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


_S3_TRIGGERING_STAGES = frozenset({"s1_rna_qc", "s2_atac_qc", "s3_doublets"})
_RNA_BRANCHES = frozenset({"paired", "unpaired", "rna_only"})


def qc_downstream_targets(run_dir: Path, stage: str) -> list[Path]:
    """The artifacts a revise of `stage` invalidates (whether or not they exist).

    The revised QC stage's own outputs, every strictly-downstream QC stage's
    outputs, AND the `post_qc_review_propose` gate outputs (proposal,
    awaiting_approval sentinel, qc_review_<run>.md, qc_summary_<run>.html). The
    gate outputs are always included: while they exist Snakemake reports "Nothing
    to be done" and silently skips the re-run even though upstream was invalidated.

    Also clears `post_qc_review.approved` when it exists — in the reprocess case
    (gate was previously approved) the sentinel must be deleted or Snakemake skips
    the QC review pause entirely and runs straight to S4-S8.

    When the revised stage triggers S3 and a post-cleanup reprocess is detected
    (rna_qc.h5ad absent on an RNA-carrying branch), the S1/S1a durable markers are
    added so Snakemake re-runs those stages rather than assuming they are done.

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
        rp.approved_sentinel("post_qc_review"),
        rp.qc_review_summary_md,
        rp.qc_summary_html,
        rp.post_qc_h5mu,
        rp.post_qc_manifest_json,
        rp.post_qc_peaks_bed,
    ]
    if stage in _S3_TRIGGERING_STAGES:
        branch = provenance.current_branch(rp.parameters_yaml)
        if branch in _RNA_BRANCHES:
            if not rp.artifact("s1_rna_qc", "rna_qc.h5ad").exists():
                targets.append(rp.artifact("s1_rna_qc", "qc_summary.json"))
            if not rp.artifact("s1a_ambient", "rna_decontaminated.h5ad").exists():
                targets.append(rp.artifact("s1a_ambient", "summary.json"))
    return targets


def invalidate_qc_downstream(run_dir: Path, stage: str) -> list[str]:
    """Delete the stale artifacts from `qc_downstream_targets` so a re-run from
    `stage` actually re-executes. Deleting non-existent files is a no-op (safe to
    call before QC has ever run, e.g. a marker-gene revise at plan review).
    Returns the paths actually deleted (for transparent logging).
    """
    deleted: list[str] = []
    for p in qc_downstream_targets(run_dir, stage):
        if p.exists():
            p.unlink()
            deleted.append(str(p))
    return deleted


def echo_binding_constraint(paths: RunPaths, stage: str) -> None:
    """Print the current QC thresholds for a QC stage so the user can see which
    bound (MAD-derived vs floor vs ceiling) is actually binding before revising —
    the pct_mt floor-vs-ceiling case is the classic gotcha."""
    if stage not in ("s1_rna_qc", "s2_atac_qc"):
        return
    import json as _json
    # Prefer LIVE effective thresholds (current parameters.yaml overlay) over the
    # frozen qc_explore.json snapshot — the snapshot is only refreshed at S0 /
    # plan-review, so it lags a post-QC revise (e.g. would show the old pct_mt
    # ceiling). Fall back to the snapshot only when the metrics parquet is absent.
    th = None
    try:
        from executor import qc_explore as _qe
        th = _qe.effective_thresholds(paths.run_dir, stage)
    except Exception:
        th = None
    if not th:
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


def revise_dry_run(run_dir: Path, paths: RunPaths, stage: str, key: str, value_parsed) -> None:
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
        existing = [p for p in qc_downstream_targets(run_dir, stage) if p.exists()]
        if existing:
            click.echo(f"  checkpoint: post_qc_review -> WOULD DELETE {len(existing)} artifact(s):")
            for p in existing:
                click.echo(f"    - {p}")
            click.echo("  Re-running the affected stages regenerates them. Confirm before "
                       "running the real `revise`.")
        else:
            click.echo("  checkpoint: post_qc_review -> no existing downstream artifacts to delete.")
    echo_binding_constraint(paths, stage)


def regenerate_plan_deliverables(run_dir: Path) -> list[str]:
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
