"""Gate enforcement + approval seeding — the mandatory pre-compute checks and the
auto-approve plumbing shared by `run` and `submit`.

Three gates fire before any compute launches: the biological-context gate, the
user-confirmed execution-mode gate, and the marker-gene decision gate. The approval
seeding helpers pre-stamp checkpoint sentinels for unattended (`--auto-approve`) batches
and resolve the Snakemake target for a resume submission.
"""
from __future__ import annotations

import os
from pathlib import Path

import click

from . import approval, context as _ctx, plan_review as _pr, provenance, stage_progress as _sp
from .cleanup import cleanup_qc_intermediates
from .log import log_event
from .pipeline import display_stage as _display_stage
from .run_paths import RunPaths


MARKER_GENE_GATE_MSG = (
    "Ambient RNA correction is planned but no marker genes are set and no explicit "
    "decision was recorded. Ask the user whether to check marker-gene expression "
    "before vs after correction (recommended, especially at elevated contamination) "
    "and to provide 5-10 gene symbols. Then re-run approve with one of:\n"
    "  - provide genes:  executor revise s1a_ambient s1a_ambient.marker_genes=\"[GENE1, GENE2]\" --config <cfg>\n"
    "  - defer to QC review:  approve plan_review --defer-marker-genes\n"
    "  - decline:  approve plan_review --skip-marker-genes\n"
    "Never invent or suggest gene symbols yourself."
)


def apply_marker_gene_ack(run_dir: Path, ack: str | None) -> None:
    """Record an unattended-batch marker-gene decision (`--marker-genes defer|skip`)."""
    if not ack:
        return
    decision = {"defer": "deferred_to_qc", "skip": "declined"}[ack]
    _pr.record_marker_gene_decision(run_dir, decision)
    log_event(run_dir, {"stage": "plan_review", "event": "marker_gene_decision",
                        "decision": decision})


def resolve_marker_gene_gate(run_dir: Path, *, defer: bool, skip: bool) -> None:
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
        raise click.ClickException(MARKER_GENE_GATE_MSG)


def enforce_context_gate(paths: RunPaths, no_context: bool) -> None:
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


def enforce_execution_mode_gate(run_dir: Path, paths: RunPaths) -> None:
    """Require an explicit, user-confirmed execution mode before launching compute.

    System requirement: Processing-MuAgent must ALWAYS confirm execution mode
    (local vs HPC) with the user before running ANY compute job — not only at S0.
    This is enforced unconditionally on every `run` and `submit`, so it also
    covers resume submissions (S1+) and runs whose config never recorded an
    execution mode. Mirrors `enforce_context_gate`.

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


def infer_submit_target(run_dir: Path) -> str:
    """Pick the Snakemake target from the first incomplete pipeline step."""
    return _sp.infer_resume_target(run_dir)


def missing_approvals(run_dir: Path, stages: tuple[str, ...]) -> list[str]:
    paths = RunPaths(run_dir)
    return [stage for stage in stages if not paths.approved_sentinel(stage).exists()]


def seed_approvals(
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
                "Cannot auto-approve plan_review: " + MARKER_GENE_GATE_MSG
                + "\nFor an unattended batch, pass --marker-genes defer|skip.")
        approval.approve(run_dir, stage, note=note)
        seeded.append(stage)
        protected = True
        if stage == "post_qc_review":
            deleted = cleanup_qc_intermediates(run_dir)
            if deleted:
                log_event(run_dir, {"stage": "post_qc_review", "event": "qc_cleanup",
                                     "deleted": deleted})
    # Set the revoke-protection flag whenever a gate is in the approved state, not
    # only when one was freshly seeded — otherwise a submit that re-enters an
    # already-approved phase would let the propose rules revoke those approvals.
    if protected:
        os.environ["PMA_AUTO_APPROVE"] = "1"
    return seeded


def prepare_submit_approvals(
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
    missing_human = missing_approvals(run_dir, human)
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
        seeded = seed_approvals(
            run_dir,
            tuple(s for s in internal if s not in _sp.required_human_approvals("all")),
            note=f"phase-auto-approved for {target}",
            kept=kept,
        )
        return seeded

    missing_internal = missing_approvals(run_dir, internal)
    if missing_internal:
        raise click.ClickException(
            f"Explicit target {target!r} requires internal approval sentinel(s): "
            f"{', '.join(missing_internal)}. Use --auto-approve for an unattended "
            "batch, approve those stages, or omit --target so submit can infer and "
            "prepare the current phase."
        )
    return []
