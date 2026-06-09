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

EXECUTOR_CHOICE = click.Choice(["local", "pbs", "slurm"])

STAGES = ["p1_context", "p2_plan", "plan_review", "s0_ingest",
          "s1a_ambient", "s1_rna_qc", "s2_atac_qc", "s3_doublets",
          "post_qc_review",
          "s4_rna_norm", "s5_atac_spectral", "s6_neighbors", "s7_clustering", "s8_umap"]

HUMAN_CHECKPOINT_STAGES = ("plan_review", "post_qc_review", "s7_clustering")
STAGE_ALIASES = {
    "qc_review": "post_qc_review",
    "resolution_review": "s7_clustering",
}
STAGE_DISPLAY = {
    "post_qc_review": "qc_review",
    "s7_clustering": "resolution_review",
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
    stage = _canonical_stage(stage)
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    _snakemake(["--configfile", str(paths.run_yaml), f"{stage}_propose"],
               run_dir, executor=executor)


def _cleanup_qc_intermediates(run_dir: Path) -> list[str]:
    """Delete large QC-only h5ad objects after post_qc_review is approved.

    Removes rna_qc.h5ad, atac_qc.h5ad, and the pre-filter atac_snap.h5ad.
    Keeps qc_summary.json, qc_metrics parquets, CBF fragment caches, S1a output,
    and all S3+ artifacts so threshold revision and downstream stages are unaffected.
    """
    rp = RunPaths(run_dir)
    targets = [
        rp.artifact("s1_rna_qc",  "rna_qc.h5ad"),
        rp.artifact("s2_atac_qc", "atac_qc.h5ad"),
        rp.artifact("s2_atac_qc", "atac_snap.h5ad"),
    ]
    deleted: list[str] = []
    for p in targets:
        if p.exists():
            p.unlink()
            deleted.append(str(p))
    return deleted


@main.command()
@click.argument("stage")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--note", default="")
def approve(stage: str, config_path: str, note: str) -> None:
    """Write internal/checkpoints/<stage>.approved to unblock <stage>_execute."""
    stage = _canonical_stage(stage)
    run_dir = _resolve_run_dir(config_path)
    approval.approve(run_dir, stage, note=note)
    log_event(run_dir, {"stage": stage, "event": "approved", "note": note})
    if stage == "post_qc_review":
        deleted = _cleanup_qc_intermediates(run_dir)
        if deleted:
            log_event(run_dir, {"stage": "post_qc_review", "event": "qc_cleanup",
                                 "deleted": deleted})
            click.echo(f"Cleaned up {len(deleted)} intermediate QC object(s).")
    click.echo(f"Approved {_display_stage(stage)}")


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
    resources_scale: float | None,
    conda_env: str | None,
) -> None:
    """Record execution mode and write deliverables/pre_run/config/site.config + hpc.env."""
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
        "resources_scale": (
            str(int(resources_scale)) if resources_scale is not None
            else os.environ.get("PMA_RESOURCES_SCALE")
        ),
        "conda_env": conda_env or os.environ.get("PMA_CONDA_ENV") or os.environ.get("CONDA_DEFAULT_ENV"),
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

    site_cfg = hpc.write_site_config(paths.site_config, mode=mode, settings=settings)
    out = hpc.write_hpc_env(paths.hpc_env_sh, paths.site_config)
    log_event(run_dir, {"stage": "configure_execution", "event": "configured",
                        "mode": mode, "hpc_env": str(out), "site_config": str(site_cfg)})
    click.echo(f"Execution mode: {mode}")
    click.echo(f"Wrote {site_cfg}")
    click.echo(f"Wrote {out}  (derived from site.config)")
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
    stage = _canonical_stage(stage)
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
    click.echo(f"Revised {key} = {value_parsed!r}; {_display_stage(stage)} is awaiting_approval.")


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
                       "`Processing-MuAgent approve <stage>` (e.g. qc_review, resolution_review).")
            return
        if paths.run_manifest_json.exists():
            click.echo("\n→ run_manifest.json present; pipeline complete.")
            return
        time.sleep(max(2.0, interval))


@main.command(name="hpc-status")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--watch", is_flag=True,
              help="Poll until a review gate needs approval, a step fails or is cancelled, "
                   "or the pipeline completes.")
@click.option("--interval", type=float, default=15.0,
              help="Poll interval in seconds when --watch is set.")
def hpc_status(config_path: str, watch: bool, interval: float) -> None:
    """Show HPC job health (monitor state) and per-step pipeline state.

    Reads Execution-MuAgent's latest_snapshot.json for health/silence/tolerance
    and latest_submission.json for the active job, then prints per-step state.
    Use --watch to poll until something actionable happens.
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)

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
        submission = _sp.load_latest_hpc_submission(paths)
        snapshot = _sp.load_hpc_monitor_state(paths)
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

    def _print(states: list[tuple[str, str, str]]) -> None:
        click.echo("")
        click.echo("--- HPC monitor ---")
        _print_hpc_header()
        click.echo("")
        click.echo("--- Pipeline state ---")
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
            click.echo(f"\n=== {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
            _print(states)
            last = states
        else:
            snapshot = _sp.load_hpc_monitor_state(paths)
            ms = (snapshot.get("monitor_state") or {}) if snapshot else {}
            health = ms.get("health", "unknown")
            silence = ms.get("silence_intervals", "?")
            tolerance = ms.get("tolerance_n", "?")
            active = next((task for _, task, st in states if st == "in_progress"), "idle")
            click.echo(
                f"[{time.strftime('%H:%M:%S')}] {active} | {health} | silence {silence}/{tolerance}"
            )
        if any(st == "failed" for _, _, st in states):
            click.echo("\n→ a step failed; inspect logs under "
                       f"{paths.snakemake_workdir}/.snakemake/slurm_logs/ "
                       "then fix and `submit` again (resume target is inferred).")
            return
        if any(st == "cancelled" for _, _, st in states):
            click.echo("\n→ a step was cancelled by the HPC monitor; see "
                       f"{paths.run_dir / 'internal' / 'hpc_monitor' / 'latest_report.md'} "
                       "for the confirmed-dead reason. Re-`submit` to resume.")
            return
        if any(st == "awaiting_approval" for _, _, st in states):
            click.echo("\n→ a review gate is awaiting approval; review deliverables and run "
                       "`Processing-MuAgent approve <stage>` (e.g. qc_review, resolution_review).")
            return
        if paths.run_manifest_json.exists():
            click.echo("\n→ run_manifest.json present; pipeline complete.")
            return
        time.sleep(max(2.0, interval))


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
    click.echo(f"Monitor: conda run --no-capture-output -n grn python -m executor.cli hpc-status --watch --config {paths.run_yaml}")


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
    """
    run_dir = _resolve_run_dir(config_path)
    if intro_context_only:
        import json as _json
        click.echo(_json.dumps(_pr.build_intro_context(run_dir), indent=2))
        return
    text = _pr.render_merged_markdown(run_dir, intro=intro_text)
    click.echo(text)
    out = _pr.write_summary(run_dir, intro=intro_text)
    click.echo(f"\nWritten: {out}")
    # Write per-stage specs; read workflow_branch from plan if available.
    try:
        import json
        plan_path = RunPaths(run_dir).artifact("p2_plan", "preprocessing_plan.json")
        branch = "paired"
        if plan_path.exists():
            branch = json.loads(plan_path.read_text()).get("workflow_branch", "paired")
        written = _specs.write_stage_specs(run_dir, branch)
        if written:
            click.echo(f"Wrote {len(written)} stage metadata file(s) to {RunPaths(run_dir).stage_meta_dir}/")
    except Exception:
        pass  # spec writing is best-effort; never block plan-review


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
@click.argument("genes", nargs=-1, required=True)
def marker_gene_check_cmd(
    config_path: str,
    force_tsne: bool,
    genes: tuple[str, ...],
) -> None:
    """Generate before/after marker gene expression plots.

    GENES is one or more gene symbols, e.g. CD3E CD20 EPCAM.

    Loads the ambient-corrected AnnData, uses a cached t-SNE embedding when the
    cell set is unchanged (or recomputes and caches it otherwise), and produces a
    side-by-side before/after expression figure.  Run
    `Processing-MuAgent propose post_qc_review` afterwards to embed the figure
    in the QC report.
    """
    from .stages import s1a_ambient as _s1a
    run_dir = _resolve_run_dir(config_path)

    if not genes:
        raise click.UsageError("Provide at least one gene symbol.")

    gene_list = list(genes)
    click.echo(f"Checking marker genes: {', '.join(gene_list)}")
    result = _s1a.run_marker_gene_check(run_dir, gene_list, force_tsne=force_tsne)
    if result["found"]:
        click.echo(f"Plotted: {', '.join(result['found'])}")
    if result["missing"]:
        click.echo(f"Not found in data: {', '.join(result['missing'])}")
    click.echo(
        "Figure written. Run `executor propose post_qc_review --config $CFG` "
        "to refresh QC reports."
    )


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

    auto_except = tuple(_canonical_stage(s) for s in auto_except)
    if auto_approve:
        # Pre-seed approval sentinels so snakemake can run the DAG end-to-end in a
        # single invocation; --auto-approve-except keeps the listed stages gated.
        kept = set(auto_except)
        _seed_approvals(run_dir, HUMAN_CHECKPOINT_STAGES, note="auto-approved", kept=kept)
        if kept:
            click.echo(f"Auto-approved all stages except: "
                       f"{sorted(_display_stage(s) for s in kept)}. "
                       "Snakemake will stop at those gates.")
    _snakemake(["--configfile", str(paths.run_yaml), target],
               run_dir, executor=executor)


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
    for stage in stages:
        if stage in kept:
            continue
        approval.approve(run_dir, stage, note=note)
        seeded.append(stage)
        if stage == "post_qc_review":
            deleted = _cleanup_qc_intermediates(run_dir)
            if deleted:
                log_event(run_dir, {"stage": "post_qc_review", "event": "qc_cleanup",
                                     "deleted": deleted})
    if seeded:
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
@click.option("--executor", type=EXECUTOR_CHOICE, required=True,
              help="Scheduler to submit the head-job to (pbs or slurm). "
                   "Use --executor local with `run` instead for foreground runs.")
@click.option("--target", default=None,
              help="Override the Snakemake target. Omit to auto-infer the first "
                   "incomplete step (e.g. s2_atac_qc_execute, post_qc_review_propose, all).")
@click.option("--auto-approve", is_flag=True,
              help="Pre-seed all checkpoint sentinels; head-job runs unattended end-to-end.")
@click.option("--auto-approve-except", "auto_except", multiple=True,
              help="With --auto-approve, keep these gates honoured. Repeatable.")
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
def submit(config_path: str, executor: str, target: str | None,
           auto_approve: bool, auto_except: tuple[str, ...],
           output_log: str | None, unlock_stale_locks: bool,
           watch: bool) -> None:
    """Submit the snakemake runner as a scheduler head-job (PBS or SLURM).

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
                --auto-approve --auto-approve-except post_qc_review \\
                --auto-approve-except s7_clustering

        # After QC review, approve and resume (target auto-inferred):
        Processing-MuAgent approve post_qc_review --config $CFG
        Processing-MuAgent submit --config $CFG --executor slurm

        # After resolution review, approve and finish:
        Processing-MuAgent approve s7_clustering --config $CFG
        Processing-MuAgent submit --config $CFG --executor slurm
    """
    if executor == "local":
        raise click.UsageError("--executor local is for `run`, not `submit`. "
                               "Use pbs or slurm here.")
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)

    auto_except = tuple(_canonical_stage(s) for s in auto_except)
    if auto_approve:
        kept = set(auto_except)
        _seed_approvals(run_dir, HUMAN_CHECKPOINT_STAGES, note="auto-approved (submit)", kept=kept)
        if kept:
            click.echo(f"Auto-approved all stages except: "
                       f"{sorted(_display_stage(s) for s in kept)}.")
        # Tell the head-job's propose rules not to revoke pre-seeded approvals.
        os.environ["PMA_AUTO_APPROVE"] = "1"

    inferred_target = target is None
    resolved_target = target if target is not None else _infer_submit_target(run_dir)
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

    out_path = Path(output_log) if output_log else hpc.head_job_log_path(executor)

    if not paths.site_config.exists():
        raise click.ClickException(
            f"site.config not found at {paths.site_config}. "
            "Run `Processing-MuAgent configure-execution --mode slurm|pbs ...` first."
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
        job_id = ea_result.get("job_id") or hpc.last_manifest_job_id(run_dir)
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
            "\nThe supervisor daemon runs kill-on-hang protection for this job. "
            "It survives SSH disconnect (unless the site uses KillUserProcesses=yes).\n"
            f"Monitor: conda run --no-capture-output -n grn python -m executor.cli hpc-status --watch --config {paths.run_yaml}"
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
        "--rerun-incomplete", *targets, *rest,
    ]
    if executor == "local":
        cmd += ["--cores", "1"]
    else:
        profile = hpc.profile_path(executor)
        cmd += ["--profile", str(profile), "--jobs", "8"]
        cmd += hpc.snakemake_cluster_cli_args()
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
