"""processing-muagent CLI — wraps Snakemake for the interactive checkpointed workflow."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import click
import yaml

from . import approval, context as _ctx, hpc, plan_review as _pr, provenance
from .log import log_event
from .run_paths import RunPaths


PACKAGE_DIR = Path(__file__).resolve().parent.parent  # processing-MuAgent/
SNAKEFILE = PACKAGE_DIR / "workflow" / "Snakefile"

EXECUTOR_CHOICE = click.Choice(["local", "pbs", "slurm"])

STAGES = ["p1_context", "p2_plan", "plan_review", "s0_ingest",
          "s1a_ambient", "s1_rna_qc", "s2_atac_qc", "s3_doublets",
          "post_qc_review",
          "s4_rna_norm", "s5_atac_lsi", "s6_dimred", "s7_clustering", "s8_umap"]


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
@click.option("--executor", type=EXECUTOR_CHOICE, default="local",
              help="Execution backend: local (default), pbs, or slurm.")
def propose(stage: str, config_path: str, executor: str) -> None:
    """Run the <stage>_propose rule."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    _snakemake(["--configfile", str(paths.run_yaml), f"{stage}_propose"],
               run_dir, executor=executor)


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


@main.command(name="hpc-info")
def hpc_info() -> None:
    """Probe the login node for scheduler queues/partitions and current PMA_* env."""
    import json
    info = hpc.discover_site()
    click.echo(json.dumps(info, indent=2, sort_keys=True))


@main.command(name="configure-execution")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--mode", "mode", required=True, type=EXECUTOR_CHOICE,
              help="Execution backend: local, pbs, or slurm.")
@click.option("--pbs-queue", default=None, help="PBS queue name (PMA_PBS_QUEUE).")
@click.option("--pbs-project", default=None, help="PBS project code (PMA_PBS_PROJECT).")
@click.option("--slurm-partition", default=None, help="SLURM partition (PMA_SLURM_PARTITION).")
@click.option("--slurm-account", default=None, help="SLURM account (PMA_SLURM_ACCOUNT).")
@click.option("--notify-email", default=None, help="Email for batch completion (PMA_NOTIFY_EMAIL).")
@click.option("--resources-scale", default=None, type=float,
              help="Memory/walltime scale factor (PMA_RESOURCES_SCALE).")
@click.option("--conda-env", default=None, help="Conda env name for cluster jobs (PMA_CONDA_ENV).")
def configure_execution(
    config_path: str,
    mode: str,
    pbs_queue: str | None,
    pbs_project: str | None,
    slurm_partition: str | None,
    slurm_account: str | None,
    notify_email: str | None,
    resources_scale: float | None,
    conda_env: str | None,
) -> None:
    """Record execution mode and write deliverables/pre_run/config/hpc.env."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    paths.ensure()

    provenance.set_param(
        str(paths.parameters_yaml),
        "execution.mode", mode,
        source="user", confidence="high",
        rationale=f"Execution backend set via configure-execution --mode {mode}.",
    )

    settings: dict[str, str | None] = {
        "pbs_queue": pbs_queue or os.environ.get("PMA_PBS_QUEUE"),
        "pbs_project": pbs_project or os.environ.get("PMA_PBS_PROJECT"),
        "slurm_partition": slurm_partition or os.environ.get("PMA_SLURM_PARTITION"),
        "slurm_account": slurm_account or os.environ.get("PMA_SLURM_ACCOUNT"),
        "notify_email": notify_email or os.environ.get("PMA_NOTIFY_EMAIL"),
        "resources_scale": (
            str(int(resources_scale)) if resources_scale is not None
            else os.environ.get("PMA_RESOURCES_SCALE")
        ),
        "conda_env": conda_env or os.environ.get("PMA_CONDA_ENV"),
    }

    if mode == "local":
        click.echo("Execution mode: local (no hpc.env written).")
        return

    if mode == "pbs" and not settings["pbs_queue"]:
        raise click.ClickException(
            "PBS mode requires --pbs-queue or PMA_PBS_QUEUE in the environment.")
    if mode == "slurm" and not settings["slurm_partition"]:
        raise click.ClickException(
            "SLURM mode requires --slurm-partition or PMA_SLURM_PARTITION in the environment.")

    out = hpc.write_hpc_env(paths.hpc_env_sh, mode=mode, settings=settings)
    log_event(run_dir, {"stage": "configure_execution", "event": "configured",
                        "mode": mode, "hpc_env": str(out)})
    click.echo(f"Execution mode: {mode}")
    click.echo(f"Wrote {out}")
    click.echo("Source this file in your shell before submit/run on the cluster:")
    click.echo(f"  source {out}")


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


def _stage_states(paths: RunPaths) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for s in STAGES:
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
        out.append((s, state))
    return out


@main.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--watch", is_flag=True,
              help="Poll until a checkpoint needs approval or manifest completes.")
@click.option("--interval", type=float, default=15.0,
              help="Poll interval in seconds when --watch is set.")
def status(config_path: str, watch: bool, interval: float) -> None:
    """Print per-stage state. With --watch, polls until something changes."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)

    def _print(states: list[tuple[str, str]]) -> None:
        for s, st in states:
            click.echo(f"  {s:20s}  {st}")

    if not watch:
        _print(_stage_states(paths))
        return

    last: list[tuple[str, str]] | None = None
    while True:
        states = _stage_states(paths)
        if states != last:
            click.echo(f"--- {time.strftime('%Y-%m-%d %H:%M:%S')} ---")
            _print(states)
            last = states
        # Stop if any stage is awaiting approval (user must act) or manifest exists.
        if any(st == "awaiting_approval" for _, st in states):
            click.echo("\n→ a stage is awaiting approval; review and run "
                       "`processing-muagent approve <stage>`.")
            return
        if paths.run_manifest_json.exists():
            click.echo("\n→ run_manifest.json present; pipeline complete.")
            return
        time.sleep(max(2.0, interval))


@main.command(name="plan-review")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def plan_review_cmd(config_path: str) -> None:
    """Render and write the merged plan-review markdown (summary + appendix)."""
    run_dir = _resolve_run_dir(config_path)
    text = _pr.render_merged_markdown(run_dir)
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
    # Refresh deliverables layout manifest after writing comparison figures.
    if (run_dir / "deliverables").exists():
        _layout.reorganise(run_dir)
    click.echo(f"Comparison written: {out}")
    click.echo(out.read_text())


@main.command(name="run")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--auto-approve", is_flag=True, help="Auto-approve every checkpoint (noninteractive).")
@click.option("--auto-approve-except", "auto_except", multiple=True,
              help="With --auto-approve, do NOT pre-seed the given stage(s). Repeatable. "
                   "Example: --auto-approve-except s7_clustering")
@click.option("--no-context", is_flag=True, help="Explicit user choice to proceed without biological context; fields marked status=missing.")
@click.option("--target", default="all")
@click.option("--executor", type=EXECUTOR_CHOICE, default="local",
              help="Execution backend: local (default), pbs, or slurm. "
                   "When pbs/slurm, snakemake stays in the foreground on this host "
                   "and dispatches per-rule cluster jobs.")
def run_pipeline(config_path: str, auto_approve: bool, auto_except: tuple[str, ...],
                 no_context: bool, target: str, executor: str) -> None:
    """Run the full DAG. With --auto-approve, checkpoints are unblocked automatically.

    Use --auto-approve-except <stage> to keep specific gates honoured (e.g. the
    S7 clustering-resolution review in headless HPC mode).
    """
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
        # Pre-seed approval sentinels so snakemake can run the DAG end-to-end in a
        # single invocation; --auto-approve-except keeps the listed stages gated.
        kept = set(auto_except)
        for s in STAGES:
            if s in kept:
                continue
            approval.approve(run_dir, s, note="auto-approved")
        os.environ["PMA_AUTO_APPROVE"] = "1"
        if kept:
            click.echo(f"Auto-approved all stages except: {sorted(kept)}. "
                       "Snakemake will stop at those gates.")
    _snakemake(["--configfile", str(paths.run_yaml), target],
               run_dir, executor=executor)


@main.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--executor", type=EXECUTOR_CHOICE, required=True,
              help="Scheduler to submit the head-job to (pbs or slurm). "
                   "Use --executor local with `run` instead for foreground runs.")
@click.option("--target", default="all",
              help="Snakemake target the head-job will run (default: all).")
@click.option("--auto-approve", is_flag=True,
              help="Pre-seed all checkpoint sentinels; head-job runs unattended end-to-end.")
@click.option("--auto-approve-except", "auto_except", multiple=True,
              help="With --auto-approve, keep these gates honoured. Repeatable.")
@click.option("--output", "output_log", type=click.Path(), default=None,
              help="Scheduler output-log path for the head-job (optional).")
def submit(config_path: str, executor: str, target: str,
           auto_approve: bool, auto_except: tuple[str, ...],
           output_log: str | None) -> None:
    """Submit the snakemake runner as a scheduler head-job (PBS or SLURM).

    The head-job runs on a compute node, activates the project conda env, and
    invokes snakemake with the cluster profile. snakemake then submits per-stage
    child jobs. The head-job exits when the DAG completes or stops at a missing
    approval gate; it emails $PMA_NOTIFY_EMAIL on exit if set.

    Typical headless workflow on HPC:

        # Run planning interactively (Phase A), then submit the heavy middle:
        processing-muagent submit --config $CFG --executor pbs \\
                --auto-approve --auto-approve-except s7_clustering

        # After the email arrives, review resolution_review.html, approve, resume:
        processing-muagent approve s7_clustering --config $CFG
        processing-muagent submit --config $CFG --executor pbs
    """
    if executor == "local":
        raise click.UsageError("--executor local is for `run`, not `submit`. "
                               "Use pbs or slurm here.")
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)

    if auto_approve:
        kept = set(auto_except)
        for s in STAGES:
            if s in kept:
                continue
            approval.approve(run_dir, s, note="auto-approved (submit)")
        if kept:
            click.echo(f"Auto-approved all stages except: {sorted(kept)}.")

    out_path = Path(output_log) if output_log else None
    job_id = hpc.submit_head_job(executor, paths.run_yaml, target=target,
                                  output_log=out_path)
    log_event(run_dir, {"stage": "submit", "event": "head_job_submitted",
                        "executor": executor, "target": target,
                        "job_id": job_id, "auto_approve": auto_approve,
                        "kept_gates": sorted(set(auto_except))})
    click.echo(f"Submitted {executor} head-job: {job_id}")
    click.echo(f"  config:  {paths.run_yaml}")
    click.echo(f"  target:  {target}")
    if os.environ.get("PMA_NOTIFY_EMAIL"):
        click.echo(f"  notify:  {os.environ['PMA_NOTIFY_EMAIL']}")
    click.echo("\nPoll progress with: processing-muagent status --watch "
               f"--config {paths.run_yaml}")


def _snakemake(args: list[str], run_dir: Path, *, executor: str = "local") -> None:
    """Invoke snakemake.

    Local mode runs snakemake with --cores 1 for reproducibility. Cluster modes
    (pbs/slurm) attach the appropriate Snakemake profile so each non-local rule
    is dispatched as a scheduler job; planning/propose/manifest rules remain
    local per the `localrules:` directive in the Snakefile.

    Expected args shape from callers: ["--configfile", <path>, <target>].
    """
    env = os.environ.copy()
    env.setdefault("PYTHONPATH", str(PACKAGE_DIR))
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

    cmd = [sys.executable, "-m", "snakemake", "-s", str(SNAKEFILE),
           "--rerun-incomplete", *targets, *rest]
    if executor == "local":
        cmd += ["--cores", "1"]
    else:
        profile = hpc.profile_path(executor)
        cmd += ["--profile", str(profile), "--jobs", "8"]
        # SLURM site-specific defaults from env vars.
        if executor == "slurm":
            if env.get("PMA_SLURM_PARTITION"):
                cmd += ["--default-resources", f"slurm_partition={env['PMA_SLURM_PARTITION']}"]
            if env.get("PMA_SLURM_ACCOUNT"):
                cmd += ["--default-resources", f"slurm_account={env['PMA_SLURM_ACCOUNT']}"]

    if configfile_path:
        cmd += ["--configfile", configfile_path]
    click.echo(f"$ {' '.join(cmd)}")
    r = subprocess.run(cmd, env=env, cwd=str(PACKAGE_DIR))
    if r.returncode != 0:
        raise click.ClickException(f"snakemake exited with {r.returncode}")


if __name__ == "__main__":
    main()
