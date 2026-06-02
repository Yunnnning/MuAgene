"""Command line interface for Execution-MuAgent."""
from __future__ import annotations

import json
from pathlib import Path

import click

from .monitor import (
    MonitorFinding,
    MonitorState,
    Submission,
    append_execution_manifest,
    load_latest_submission,
    load_site_config,
    load_stage_spec,
    monitor_once,
    monitor_watch,
    parse_job_ids_from_log,
    register_submission,
    run_monitor_dir,
    submit_from_spec,
    utc_now,
    validate_spec,
)


@click.group()
def main() -> None:
    """Central execution monitor for MuAgene subagents."""


@main.command()
@click.option("--agent", required=True, help="Submitting subagent name, e.g. Processing-MuAgent.")
@click.option("--executor", required=True, type=click.Choice(["slurm", "pbs"]))
@click.option("--job-id", required=True, help="Scheduler head-job id.")
@click.option("--run-dir", required=True, type=click.Path())
@click.option("--config", required=True, type=click.Path())
@click.option("--target", required=True)
@click.option("--repo-root", required=True, type=click.Path())
@click.option("--log-path", required=True, type=click.Path())
@click.option("--spec-path", default=None, type=click.Path(),
              help="Path to the head-job spec YAML. When set, the monitor reads "
                   "progress_timeout_hint from the spec instead of --stale-minutes.")
def register(
    agent: str,
    executor: str,
    job_id: str,
    run_dir: str,
    config: str,
    target: str,
    repo_root: str,
    log_path: str,
    spec_path: str | None,
) -> None:
    """Break-glass: register a manually-submitted HPC job for monitoring.

    Use this when you submitted the head-job manually with sbatch/qsub (instead of
    via `Processing-MuAgent submit`). Execution-MuAgent will pick up the job and
    monitor it as if it had submitted the job itself.
    """
    progress_timeout: float | None = None
    if spec_path and Path(spec_path).exists():
        try:
            progress_timeout = load_stage_spec(spec_path).progress_timeout_hint
        except Exception:
            pass
    submission = Submission(
        agent=agent,
        executor=executor,
        job_id=job_id,
        run_dir=str(Path(run_dir).resolve()),
        config=str(Path(config).resolve()),
        target=target,
        repo_root=str(Path(repo_root).resolve()),
        log_path=str(Path(log_path).resolve()),
        submitted_at=utc_now(),
        spec_path=str(Path(spec_path).resolve()) if spec_path else None,
        progress_timeout_hint=progress_timeout,
    )
    path = register_submission(submission)
    click.echo(f"registered {executor} job {job_id} for {agent}")
    click.echo(f"registry: {path}")


@main.command()
@click.option("--run-dir", required=True, type=click.Path(exists=True))
@click.option("--job-id", default=None, help="Specific registered job id. Defaults to latest.")
@click.option("--stale-minutes", default=90.0, show_default=True, type=float,
              help="Fallback timeout (minutes) used to derive tolerance_n when no spec is linked. "
                   "Spec's progress_timeout_hint takes priority.")
@click.option("--scheduler-timeout", default=5, show_default=True, type=int)
@click.option("--kill-on-hang/--no-kill-on-hang", default=True, show_default=True,
              help="Cancel child jobs first, then the head job, on confirmed-dead classification.")
@click.option("--watch", is_flag=True, help="Continue polling until all jobs finish.")
@click.option("--interval", default=900.0, show_default=True, type=float,
              help="Check interval in seconds. tolerance_n = ceil(progress_timeout / interval).")
@click.option("--max-checks", default=None, type=int, help="Stop after N checks when --watch is used.")
@click.option("--json-output", is_flag=True, help="Print machine-readable snapshot.")
@click.option("--spec-path", default=None, type=click.Path(),
              help="Per-stage spec YAML; overrides --stale-minutes with spec's progress_timeout_hint.")
@click.option("--fs-hang-policy", default="hold", show_default=True,
              type=click.Choice(["hold", "kill_and_resubmit"]),
              help="What to do when a filesystem hang is detected. "
                   "'hold' reports and waits; 'kill_and_resubmit' cancels and routes through failure path.")
def monitor(
    run_dir: str,
    job_id: str | None,
    stale_minutes: float,
    scheduler_timeout: int,
    kill_on_hang: bool,
    watch: bool,
    interval: float,
    max_checks: int | None,
    json_output: bool,
    spec_path: str | None,
    fs_hang_policy: str,
) -> None:
    """Monitor a registered HPC job and write a report when problems are detected.

    The monitor runs a two-clock state machine: the watcher counts quiet check
    intervals (no file progress, no log growth); after tolerance_n quiet intervals
    the job enters SUSPECT and investigation begins. Kill only from CONFIRMED_DEAD.
    """
    submission = load_latest_submission(run_dir, job_id=job_id)
    if spec_path and not submission.spec_path:
        import dataclasses
        submission = dataclasses.replace(submission, spec_path=str(Path(spec_path).resolve()))
    if watch:
        report = monitor_watch(
            submission,
            interval_s=interval,
            stale_minutes=stale_minutes,
            scheduler_timeout_s=scheduler_timeout,
            kill_on_hang=kill_on_hang,
            fs_hang_policy=fs_hang_policy,
            max_checks=max_checks,
        )
        if report:
            click.echo(f"problem report: {report}")
        else:
            click.echo("no problem report written")
        return
    from .monitor import _resolve_tolerance_n
    tolerance_n = _resolve_tolerance_n(submission, interval, stale_minutes)
    state = MonitorState(tolerance_n=tolerance_n)
    snapshot, findings, cancel_result, report, state = monitor_once(
        submission,
        state,
        scheduler_timeout_s=scheduler_timeout,
        kill_on_hang=kill_on_hang,
        fs_hang_policy=fs_hang_policy,
    )
    if json_output:
        click.echo(json.dumps({
            "snapshot": snapshot,
            "findings": [f.__dict__ for f in findings],
            "cancel_result": cancel_result,
            "report": str(report) if report else None,
            "monitor_state": {"health": state.health.value, "silence_intervals": state.silence_intervals},
        }, indent=2, sort_keys=True))
    else:
        if findings:
            click.echo(f"{len(findings)} finding(s) detected")
            for finding in findings:
                click.echo(f"  {finding.severity}: {finding.code}: {finding.message}")
        else:
            click.echo("no monitor problems detected")
        click.echo(f"health: {state.health.value}  silence: {state.silence_intervals}/{state.tolerance_n}")
        if cancel_result:
            click.echo(f"cancel attempted: {cancel_result}")
        if report:
            click.echo(f"problem report: {report}")


@main.command(name="execute-spec")
@click.option("--spec", "spec_path", required=True, type=click.Path(exists=True),
              help="Path to the head-job spec YAML (internal/stage_meta/head_job.yaml).")
@click.option("--site-config", "site_config_path", required=True, type=click.Path(exists=True),
              help="Path to site.config YAML (platform description).")
@click.option("--run-dir", required=True, type=click.Path())
@click.option("--repo-root", required=True, type=click.Path(exists=True),
              help="Processing-MuAgent repo root (for runner scripts).")
@click.option("--target", "target_arg", default=None,
              help="Snakemake target to pass as PMA_TARGET (e.g. 'all', 'post_qc_review_propose'). "
                   "Defaults to '<spec.stage>_execute' when omitted.")
@click.option("--watch", is_flag=True, help="Monitor the job until it exits after submission.")
@click.option("--interval", default=900.0, show_default=True, type=float,
              help="Check interval in seconds when --watch is used.")
@click.option("--kill-on-hang/--no-kill-on-hang", default=True, show_default=True)
def execute_spec(
    spec_path: str,
    site_config_path: str,
    run_dir: str,
    repo_root: str,
    target_arg: str | None,
    watch: bool,
    interval: float,
    kill_on_hang: bool,
) -> None:
    """Validate the head-job spec, render a submission script, submit, record, and optionally monitor.

    Execution-MuAgent owns everything between a spec and a running job. Processing-MuAgent
    writes the head_job.yaml spec and site.config; this command handles the rest without
    user interaction. Snakemake submits per-stage child jobs from within the head-job.
    """
    spec = load_stage_spec(spec_path)
    site_cfg = load_site_config(site_config_path)

    errors = validate_spec(spec, site_cfg)
    if errors:
        for err in errors:
            click.echo(f"validation error: {err}", err=True)
        raise click.ClickException(
            f"Spec validation failed with {len(errors)} error(s); not submitting."
        )

    run_dir_path = Path(run_dir).resolve()
    log_dir = run_dir_path / "internal" / "hpc_monitor" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp = utc_now().replace(":", "").replace("-", "")
    log_path = log_dir / f"{spec.stage}_{stamp}.out"
    target = target_arg or f"{spec.stage}_execute"

    result = submit_from_spec(
        spec, site_cfg, run_dir_path, repo_root, log_path, target,
    )

    if result["rejected_as"] == "policy":
        finding = MonitorFinding(
            severity="error",
            code="submit_rejected_policy",
            message=(
                f"Submission of {spec.stage!r} was rejected by the scheduler as a policy error "
                f"(invalid partition/account/walltime). Scheduler said: {result['stderr'] or result['stdout']}. "
                "Adjust partition, account, or resources_scale in site.config and resubmit."
            ),
        )
        from .monitor import render_report, write_report, collect_snapshot
        report_text = (
            f"# HPC Monitor Report\n\n"
            f"## Findings\n\n"
            f"- **ERROR `{finding.code}`**: {finding.message}\n"
        )
        latest = run_monitor_dir(run_dir_path) / "latest_report.md"
        latest.parent.mkdir(parents=True, exist_ok=True)
        latest.write_text(report_text, encoding="utf-8")
        click.echo(f"submit_rejected_policy: {finding.message}", err=True)
        raise click.ClickException("Submission rejected (policy); see latest_report.md.")

    if result["rejected_as"] == "transient":
        raise click.ClickException(
            f"Submission of {spec.stage!r} failed after retries (transient scheduler error). "
            f"stderr: {result['stderr']}"
        )

    job_id = result["job_id"]
    click.echo(f"Submitted {site_cfg.scheduler} head-job: {job_id}")
    click.echo(f"  stage:  {spec.stage}")
    click.echo(f"  script: {result['script_path']}")
    click.echo(f"  log:    {log_path}")

    append_execution_manifest(run_dir_path, {
        "submitted_at": utc_now(),
        "stage": spec.stage,
        "science_description": spec.science_description,
        "job_id": job_id,
        "spec_path": str(Path(spec_path).resolve()),
        "script_path": result["script_path"],
        "run_dir": str(run_dir_path),
        "expected_outputs": spec.outputs,
    })

    submission = Submission(
        agent="Execution-MuAgent",
        executor=site_cfg.scheduler,
        job_id=job_id,
        run_dir=str(run_dir_path),
        config=str(run_dir_path / "deliverables" / "pre_run" / "config" / "run.yaml"),
        target=target,
        repo_root=str(Path(repo_root).resolve()),
        log_path=str(log_path),
        submitted_at=utc_now(),
        spec_path=str(Path(spec_path).resolve()),
        progress_timeout_hint=spec.progress_timeout_hint,
    )
    registry = register_submission(submission)
    click.echo(f"  monitor registry: {registry}")

    if watch:
        report = monitor_watch(
            submission,
            interval_s=interval,
            stale_minutes=spec.progress_timeout_hint,
            scheduler_timeout_s=5,
            kill_on_hang=kill_on_hang,
            fs_hang_policy=site_cfg.fs_hang_policy,
            max_checks=None,
        )
        if report:
            click.echo(f"problem report: {report}")
        else:
            click.echo("no problem report written")


@main.command()
@click.argument("log_path", type=click.Path(exists=True))
def parse_jobs(log_path: str) -> None:
    """Extract scheduler job ids from a runner or Snakemake log."""
    for job_id in parse_job_ids_from_log(log_path):
        click.echo(job_id)


if __name__ == "__main__":
    main()
