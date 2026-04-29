"""processing-muagent CLI — wraps Snakemake for the interactive checkpointed workflow."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import click
import yaml

from . import approval, context as _ctx, plan_review as _pr, provenance
from .log import log_event
from .run_paths import RunPaths


PACKAGE_DIR = Path(__file__).resolve().parent.parent  # processing-MuAgent/
SNAKEFILE = PACKAGE_DIR / "workflow" / "Snakefile"


def _resolve_run_dir(config_path: Path | str) -> Path:
    with Path(config_path).open() as f:
        cfg = yaml.safe_load(f) or {}
    rd = cfg.get("run_dir")
    if not rd:
        raise click.ClickException("run.yaml must set 'run_dir'")
    return Path(rd).expanduser().resolve()


@click.group()
def main() -> None:
    """processing-MuAgent: multiome preprocessing subagent (stops after per-modality UMAP)."""


@main.command()
@click.option("--config", "config_path", required=True, type=click.Path())
def init(config_path: str) -> None:
    """Initialize a run directory.

    Creates the `internal/` and `deliverables/` scaffolds, copies the user's
    config into its canonical user-facing location `deliverables/config/run.yaml`,
    and writes the Biological Context Report template into
    `deliverables/config/biological_context.md`.
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
    """Run the <stage>_propose rule."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    _snakemake(["--configfile", str(paths.run_yaml), f"{stage}_propose"], run_dir)


@main.command()
@click.argument("stage")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--note", default="")
def approve(stage: str, config_path: str, note: str) -> None:
    """Write internal/checkpoints/<stage>.approved to unblock <stage>_execute."""
    run_dir = _resolve_run_dir(config_path)
    approval.approve(run_dir, stage, note=note)
    log_event(run_dir, {"stage": stage, "event": "approved", "note": note})
    click.echo(f"Approved {stage}")


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


@main.command()
@click.argument("stage")
@click.argument("param_kv")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--rationale", default="User revision")
def revise(stage: str, param_kv: str, config_path: str, rationale: str) -> None:
    """Update one parameter and reset the stage to awaiting_approval.

    PARAM_KV is key=value, e.g. s1_rna_qc.pct_counts_mt_max=10.0
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    if "=" not in param_kv:
        raise click.ClickException("param_kv must be key=value")
    key, value = param_kv.split("=", 1)
    try:
        value_parsed = yaml.safe_load(value)
    except Exception:
        value_parsed = value
    provenance.set_param(
        str(paths.parameters_yaml),
        key, value_parsed,
        source="user", confidence="high", rationale=rationale,
    )
    approval.mark_awaiting(run_dir, stage)
    log_event(run_dir, {"stage": stage, "event": "revised", "param": key, "value": value_parsed})
    click.echo(f"Revised {key} = {value_parsed!r}; {stage} is awaiting_approval.")


@main.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def status(config_path: str) -> None:
    """Print per-stage state."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    stages = ["p1_context", "p2_plan", "plan_review", "s0_ingest",
              "s1a_ambient", "s1_rna_qc", "s2_atac_qc", "s3_doublets", "s4_rna_norm",
              "s5_atac_lsi", "s6_dimred", "s7_clustering", "s8_umap"]
    for s in stages:
        proposed = paths.proposal(s).exists()
        awaiting = paths.awaiting_sentinel(s).exists()
        approved = paths.approved_sentinel(s).exists()
        if approved and not awaiting:
            state = "approved"
        elif awaiting:
            state = "awaiting_approval"
        elif proposed:
            state = "proposed"
        else:
            state = "pending"
        click.echo(f"  {s:20s}  {state}")


@main.command(name="plan-review")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def plan_review_cmd(config_path: str) -> None:
    """Print the concise preprocessing-plan review summary for the run."""
    run_dir = _resolve_run_dir(config_path)
    items = _pr.build_summary(run_dir)
    text = _pr.render_summary_text(items)
    click.echo(text)
    out = _pr.write_summary(run_dir)
    click.echo(f"\nWritten: {out}")


@main.command(name="resolution-compare")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--rna", "rna_pair", default="1.0,1.2",
              help="Comma-separated RNA resolutions, e.g. 1.0,1.2")
@click.option("--atac", "atac_pair", default="0.6,0.8",
              help="Comma-separated ATAC resolutions, e.g. 0.6,0.8")
def resolution_compare_cmd(config_path: str, rna_pair: str, atac_pair: str) -> None:
    """Render side-by-side Leiden resolution comparisons for RNA and ATAC.

    Re-clusters at the specified resolutions; does NOT change the approved cluster
    labels. Produces side-by-side UMAP figures + a markdown summary.
    """
    from . import resolution_compare as _rc, layout as _layout
    run_dir = _resolve_run_dir(config_path)
    rna_res = tuple(float(x) for x in rna_pair.split(","))
    atac_res = tuple(float(x) for x in atac_pair.split(","))
    out = _rc.run_comparison(run_dir, rna_resolutions=rna_res, atac_resolutions=atac_res)
    # Refresh deliverables/figures so the new comparison PNG+PDF show up there (if they qualify).
    if (run_dir / "deliverables").exists():
        _layout.reorganise(run_dir)
    click.echo(f"Comparison written: {out}")
    click.echo(out.read_text())


@main.command(name="run")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--auto-approve", is_flag=True, help="Auto-approve every checkpoint (noninteractive).")
@click.option("--no-context", is_flag=True, help="Explicit user choice to proceed without biological context; fields marked status=missing.")
@click.option("--target", default="all")
def run_pipeline(config_path: str, auto_approve: bool, no_context: bool, target: str) -> None:
    """Run the full DAG. With --auto-approve, checkpoints are unblocked automatically."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)

    # Phase 1 biological context check (MANDATORY FIRST STEP).
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

    if auto_approve:
        # Pre-seed all approval sentinels so snakemake can run the DAG end-to-end
        # in a single invocation; set PMA_AUTO_APPROVE=1 so propose rules don't
        # strip the approvals when they (re-)run.
        stages = ["p1_context", "p2_plan", "plan_review", "s0_ingest",
                  "s1a_ambient", "s1_rna_qc", "s2_atac_qc", "s3_doublets", "s4_rna_norm",
                  "s5_atac_lsi", "s6_dimred", "s7_clustering", "s8_umap"]
        for s in stages:
            approval.approve(run_dir, s, note="auto-approved")
        os.environ["PMA_AUTO_APPROVE"] = "1"
        _snakemake(["--configfile", str(paths.run_yaml), "all"], run_dir)
    else:
        _snakemake(["--configfile", str(paths.run_yaml), target], run_dir)


def _snakemake(args: list[str], run_dir: Path) -> None:
    """Invoke snakemake. NOTE: --configfile has nargs=+, so targets must NOT appear
    directly after it. We put targets first, then options at the end.
    Expected args shape from callers: ["--configfile", <path>, <target>].
    """
    env = os.environ.copy()
    env.setdefault("PYTHONPATH", str(PACKAGE_DIR))
    # Single-thread for reproducibility (UMAP / numba).
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
    cmd = [sys.executable, "-m", "snakemake", "-s", str(SNAKEFILE),
           "--cores", "1", "--rerun-incomplete", *targets, *rest]
    if configfile_path:
        cmd += ["--configfile", configfile_path]
    click.echo(f"$ {' '.join(cmd)}")
    r = subprocess.run(cmd, env=env)
    if r.returncode != 0:
        raise click.ClickException(f"snakemake exited with {r.returncode}")


if __name__ == "__main__":
    main()
