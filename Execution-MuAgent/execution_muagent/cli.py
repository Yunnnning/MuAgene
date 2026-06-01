"""Command line interface for Execution-MuAgent."""
from __future__ import annotations

import json
from pathlib import Path

import click

from .monitor import (
    Submission,
    load_latest_submission,
    monitor_once,
    monitor_watch,
    parse_job_ids_from_log,
    register_submission,
    utc_now,
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
def register(
    agent: str,
    executor: str,
    job_id: str,
    run_dir: str,
    config: str,
    target: str,
    repo_root: str,
    log_path: str,
) -> None:
    """Register an HPC submission for monitoring."""
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
    )
    path = register_submission(submission)
    click.echo(f"registered {executor} job {job_id} for {agent}")
    click.echo(f"registry: {path}")


@main.command()
@click.option("--run-dir", required=True, type=click.Path(exists=True))
@click.option("--job-id", default=None, help="Specific registered job id. Defaults to latest.")
@click.option("--stale-minutes", default=20.0, show_default=True, type=float)
@click.option("--scheduler-timeout", default=5, show_default=True, type=int)
@click.option("--kill-on-hang", is_flag=True, help="Cancel the head job on detected hang/failure.")
@click.option("--watch", is_flag=True, help="Continue polling until a problem is reported.")
@click.option("--interval", default=60.0, show_default=True, type=float)
@click.option("--max-checks", default=None, type=int, help="Stop after N checks when --watch is used.")
@click.option("--json-output", is_flag=True, help="Print machine-readable snapshot.")
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
) -> None:
    """Monitor a registered HPC job and write a report when problems are detected."""
    submission = load_latest_submission(run_dir, job_id=job_id)
    if watch:
        report = monitor_watch(
            submission,
            interval_s=interval,
            stale_minutes=stale_minutes,
            scheduler_timeout_s=scheduler_timeout,
            kill_on_hang=kill_on_hang,
            max_checks=max_checks,
        )
        if report:
            click.echo(f"problem report: {report}")
        else:
            click.echo("no problem report written")
        return
    snapshot, findings, cancel_result, report = monitor_once(
        submission,
        stale_minutes=stale_minutes,
        scheduler_timeout_s=scheduler_timeout,
        kill_on_hang=kill_on_hang,
    )
    if json_output:
        click.echo(json.dumps({
            "snapshot": snapshot,
            "findings": [f.__dict__ for f in findings],
            "cancel_result": cancel_result,
            "report": str(report) if report else None,
        }, indent=2, sort_keys=True))
    else:
        if findings:
            click.echo(f"{len(findings)} finding(s) detected")
            for finding in findings:
                click.echo(f"  {finding.severity}: {finding.code}: {finding.message}")
        else:
            click.echo("no monitor problems detected")
        if cancel_result:
            click.echo(f"cancel attempted: {cancel_result}")
        if report:
            click.echo(f"problem report: {report}")


@main.command()
@click.argument("log_path", type=click.Path(exists=True))
def parse_jobs(log_path: str) -> None:
    """Extract scheduler job ids from a runner or Snakemake log."""
    for job_id in parse_job_ids_from_log(log_path):
        click.echo(job_id)


if __name__ == "__main__":
    main()
