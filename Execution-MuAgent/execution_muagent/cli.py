"""Command line interface for Execution-MuAgent."""
from __future__ import annotations

from pathlib import Path

import click

from .monitor import (
    MonitorFinding,
    Submission,
    append_execution_manifest,
    load_site_config,
    load_stage_spec,
    monitor_watch,
    register_submission,
    run_monitor_dir,
    submit_from_spec,
    utc_now,
    validate_spec,
)


@click.group()
def main() -> None:
    """Central execution monitor for MuAgene subagents."""


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
@click.option("--interval", default=270.0, show_default=True, type=float,
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
            max_checks=None,
        )
        if report:
            click.echo(f"problem report: {report}")
        else:
            click.echo("no problem report written")


@main.command()
@click.option("--run-dir", required=True, type=click.Path(exists=True),
              help="Run directory whose diagnostic report should be printed.")
def report(run_dir: str) -> None:
    """Print the latest diagnostic report written by Execution-MuAgent.

    The only human-facing command. Reads
    `<run_dir>/internal/hpc_monitor/latest_report.md` — the findings and
    confirmed-dead/verification diagnostics for the most recent monitor check —
    and prints it to stdout. All other behaviour (submit, monitor, verify, kill)
    is driven by Processing-MuAgent via `execute-spec`.
    """
    latest = run_monitor_dir(run_dir) / "latest_report.md"
    if not latest.is_file():
        raise click.ClickException(f"No diagnostic report at {latest}")
    click.echo(latest.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
