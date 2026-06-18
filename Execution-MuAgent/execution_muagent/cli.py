"""Command line interface for Execution-MuAgent."""
from __future__ import annotations

from pathlib import Path

import click

from .monitor import (
    DEFAULT_CHECK_INTERVAL_S,
    MonitorFinding,
    Submission,
    append_execution_manifest,
    load_site_config,
    load_stage_spec,
    monitor_watch,
    register_submission,
    resolve_config_yaml,
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
@click.option("--interval", default=DEFAULT_CHECK_INTERVAL_S, show_default=True, type=float,
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

    # Environment preflight: make the run's env(s) real + valid before submitting.
    # CPU env always (head job + CPU stages). GPU env too, but ONLY when device=gpu AND
    # this run actually has a GPU-capable stage (spec.gpu_stages_present, set by Processing
    # from _GPU_CAPABLE) — otherwise a device=gpu preprocessing run (no GPU consumer) would
    # pull the multi-GB container for nothing. Deliberate eager prep stays available via an
    # explicit `provision-env --device gpu`. policy=auto provisions a missing/stale env;
    # never silently degrade — a missing/invalid GPU env aborts the submit (it would fail
    # mid-run otherwise). Old site.configs without an `environments:` section reconcile to
    # "ok" without action, so existing CPU runs are unaffected.
    from . import environment
    gpu_needed = (site_cfg.device or "cpu") == "gpu" and getattr(spec, "gpu_stages_present", False)
    devices = ["cpu"] + (["gpu"] if gpu_needed else [])
    env_errors: list[str] = []
    for dev in devices:
        try:
            rec = environment.reconcile(site_cfg, repo_root, dev)
        except Exception as exc:
            # A crashed reconcile (corrupt env_state.json, probe/subprocess failure)
            # must NOT be downgraded to a warning — that would submit a job against an
            # unverified env (silent degrade). Record it as a hard preflight error so
            # the `if env_errors` check below aborts before submit_from_spec.
            env_errors.append(f"{dev}: env preflight crashed: {exc}")
            click.echo(f"env preflight [{dev}] error: {exc}", err=True)
            continue
        if rec.get("provision"):
            click.echo(f"env preflight [{dev}]: "
                       f"{rec['provision'].get('action')} -> {rec['provision'].get('status')}")
        for f in rec.get("findings", []):
            click.echo(f"env preflight [{dev}] {f['severity']}: {f['message']}", err=True)
            if f["severity"] == "error":
                env_errors.append(f"{dev}: {f['message']}")
    if env_errors:
        raise click.ClickException(
            "Environment preflight failed; not submitting:\n  - " + "\n  - ".join(env_errors))

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
        config=str(resolve_config_yaml(run_dir_path)),
        target=target,
        repo_root=str(Path(repo_root).resolve()),
        log_path=str(log_path),
        submitted_at=utc_now(),
        spec_path=str(Path(spec_path).resolve()),
        progress_timeout_hint=spec.progress_timeout_hint,
        device=site_cfg.device or "cpu",
    )
    registry = register_submission(submission)
    click.echo(f"  monitor registry: {registry}")

    if watch:
        pid_path = run_dir_path / "internal" / "hpc_monitor" / "monitor.pid"
        try:
            report = monitor_watch(
                submission,
                interval_s=interval,
                stale_minutes=spec.progress_timeout_hint,
                scheduler_timeout_s=5,
                kill_on_hang=kill_on_hang,
                max_checks=None,
            )
        finally:
            pid_path.unlink(missing_ok=True)
        if report:
            click.echo(f"problem report: {report}")
        else:
            click.echo("no problem report written")


@main.command(name="resume-monitor")
@click.option("--run-dir", required=True, type=click.Path())
@click.option("--interval", default=DEFAULT_CHECK_INTERVAL_S, show_default=True, type=float,
              help="Check interval in seconds.")
@click.option("--kill-on-hang/--no-kill-on-hang", default=True, show_default=True)
def resume_monitor(run_dir: str, interval: float, kill_on_hang: bool) -> None:
    """Resume monitoring an existing submission without resubmitting the cluster job.

    Reads latest_submission.json from the run's hpc_monitor directory, reconstructs
    the Submission dataclass, and calls monitor_watch(). Intended to be invoked as a
    background daemon by Processing-MuAgent supervisor-restart. All output is captured
    by the caller (goes to monitor_<ts>.log — never the user terminal).
    """
    import dataclasses
    import json

    run_dir_path = Path(run_dir).resolve()
    sub_path = run_dir_path / "internal" / "hpc_monitor" / "latest_submission.json"
    if not sub_path.exists():
        raise click.ClickException(f"No submission found: {sub_path}")
    data = json.loads(sub_path.read_text())
    field_names = {f.name for f in dataclasses.fields(Submission)}
    sub = Submission(**{k: v for k, v in data.items() if k in field_names})
    stale_minutes = sub.progress_timeout_hint or 90.0
    pid_path = run_dir_path / "internal" / "hpc_monitor" / "monitor.pid"
    try:
        report_path = monitor_watch(
            sub,
            interval_s=interval,
            stale_minutes=stale_minutes,
            scheduler_timeout_s=5,
            kill_on_hang=kill_on_hang,
        )
    finally:
        pid_path.unlink(missing_ok=True)
    if report_path:
        click.echo(report_path.read_text())
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


def _provisioning_context(site_config_path: str | None, repo_root: str | None):
    """Resolve (site_cfg, repo_root) for provision/validate from EITHER an explicit
    science site.config OR this machine's ~/.muagene/machine.config profile.

    With --site-config: the backward-compatible path (its `environments:` section
    wins; pre-contract configs still work). Without it: synthesize a site.config from
    machine.config + the committed env manifest, so a fresh machine can provision
    before any run/site.config exists — Execution-MuAgent owns infra end to end.
    """
    from . import machine
    mc = machine.load_machine_config()
    if site_config_path:
        site_cfg = load_site_config(site_config_path)
        repo = repo_root or (mc.processing_repo if mc else None)
        if not repo:
            raise click.ClickException(
                "--repo-root is required with --site-config (the Processing-MuAgent repo "
                "holding env definitions), and no machine.config processing_repo is recorded.")
        return site_cfg, str(repo)
    repo = repo_root or (mc.processing_repo if mc else None)
    if not repo:
        raise click.ClickException(
            "No --site-config given and this machine is not bootstrapped. Run "
            "`Execution-MuAgent init-machine --processing-repo <path>` first, or pass "
            "--site-config and --repo-root explicitly.")
    return machine.synthesize_site_config(repo, mc), str(repo)


@main.command(name="provision-env")
@click.option("--site-config", "site_config_path", default=None, type=click.Path(exists=True),
              help="Optional site.config (its environments: section drives provisioning). "
                   "Omit to use this machine's ~/.muagene/machine.config profile.")
@click.option("--repo-root", default=None, type=click.Path(exists=True),
              help="Processing-MuAgent repo root (holds env definitions/locks/.def). "
                   "Defaults to machine.config processing_repo.")
@click.option("--device", type=click.Choice(["cpu", "gpu", "both"]), default="both", show_default=True)
@click.option("--force", is_flag=True, help="Re-provision even if the env is already present and current.")
def provision_env_cmd(site_config_path: str | None, repo_root: str | None, device: str, force: bool) -> None:
    """Make the CPU/GPU env real on THIS machine (idempotent).

    CPU = conda-lock create; GPU = pull a pinned, centrally-published container image
    (never built locally). Records a fingerprint so a later definition/image change is
    detected and re-provisioned. Run once per new machine (or via `init-machine`).
    """
    from . import environment
    site_cfg, repo_root = _provisioning_context(site_config_path, repo_root)
    devices = ["cpu", "gpu"] if device == "both" else [device]
    failed = False
    for dev in devices:
        spec = environment.resolve_env_spec(site_cfg, repo_root, dev)
        click.echo(f"[{dev}] provider={spec.provider} "
                   f"target={spec.image or spec.env_name}")
        res = environment.provision_env(spec, site_cfg, force=force)
        click.echo(f"[{dev}] {res.get('action', 'noop')} -> {res.get('status')}")
        if res.get("status") == "failed":
            click.echo((res.get("stderr") or "")[-1200:], err=True)
            failed = True
    if failed:
        raise click.ClickException("provision-env failed for one or more devices (see above).")


@main.command(name="validate-env")
@click.option("--site-config", "site_config_path", default=None, type=click.Path(exists=True),
              help="Optional site.config. Omit to use this machine's machine.config profile.")
@click.option("--repo-root", default=None, type=click.Path(exists=True),
              help="Processing-MuAgent repo root. Defaults to machine.config processing_repo.")
@click.option("--device", type=click.Choice(["cpu", "gpu", "both"]), default="both", show_default=True)
def validate_env_cmd(site_config_path: str | None, repo_root: str | None, device: str) -> None:
    """Preflight an env without submitting: present + imports its declared modules.

    Nonzero exit on any error finding (e.g. a CPU-only env asked to run GPU work).
    """
    from . import environment
    site_cfg, repo_root = _provisioning_context(site_config_path, repo_root)
    devices = ["cpu", "gpu"] if device == "both" else [device]
    ok = True
    for dev in devices:
        spec = environment.resolve_env_spec(site_cfg, repo_root, dev)
        res = environment.validate_env(spec, site_cfg)
        if not res["findings"]:
            click.echo(f"[{dev}] OK ({spec.provider})")
        for f in res["findings"]:
            click.echo(f"[{dev}] {f['severity']}: {f['message']}")
        ok = ok and res["ok"]
    if not ok:
        raise click.ClickException("validate-env found errors (see above).")


@main.command()
@click.option("--site-config", "site_config_path", default=None, type=click.Path(exists=True),
              help="Optional: validate the envs declared in this site.config.")
@click.option("--repo-root", default=None, type=click.Path(exists=True))
def doctor(site_config_path: str | None, repo_root: str | None) -> None:
    """Print this machine's capabilities (scheduler, GPU, env manager, container
    runtime) and validate its envs (from the given site.config, else machine.config)."""
    import json as _json

    from . import capabilities, environment, machine
    caps = capabilities.probe_capabilities()
    click.echo(_json.dumps(caps, indent=2, sort_keys=True))
    site_cfg = None
    if site_config_path:
        site_cfg = load_site_config(site_config_path)
        repo = repo_root or (machine.load_machine_config() or machine.MachineConfig()).processing_repo
    else:
        mc = machine.load_machine_config()
        repo = repo_root or (mc.processing_repo if mc else None)
        if mc and repo:
            site_cfg = machine.synthesize_site_config(repo, mc)
            click.echo(f"machine profile: {machine.machine_config_path()}")
    if site_cfg is None or not repo:
        return
    for dev in ("cpu", "gpu"):
        spec = environment.resolve_env_spec(site_cfg, repo, dev)
        res = environment.validate_env(spec, site_cfg)
        status = "OK" if res["ok"] else "PROBLEM"
        click.echo(f"[{dev}] {status} (provider={spec.provider})")
        for f in res["findings"]:
            click.echo(f"  {f['severity']}: {f['message']}")


@main.command(name="init-machine")
@click.option("--processing-repo", required=True, type=click.Path(exists=True),
              help="Sibling Processing-MuAgent repo root (holds env definitions + manifest.yaml).")
@click.option("--device", type=click.Choice(["cpu", "gpu", "both"]), default="cpu", show_default=True)
@click.option("--manager", default=None, help="conda env manager (default: auto-detect).")
@click.option("--container-runtime", default=None, help="apptainer|singularity (default: auto-detect).")
@click.option("--singularity-module", default=None, help="`module load` name for singularity on HPC.")
@click.option("--gpu-image", default=None, help="Machine-local .sif path the image pulls to.")
@click.option("--gpu-image-uri", default=None,
              help="Pinned registry reference to PULL the GPU image from (required when device "
                   "includes gpu). No machine builds the image locally.")
@click.option("--conda-env", default=None, help="CPU env name to create (default muagene).")
@click.option("--gpu-conda-env", default=None, help="GPU env name (default muagene-gpu).")
@click.option("--policy", type=click.Choice(["auto", "manual"]), default="auto", show_default=True,
              help="Submit-time reconcile policy recorded in machine.config.")
@click.option("--install-processing/--no-install-processing", default=True, show_default=True,
              help="pip install -e the Processing + Execution packages into the CPU env.")
@click.option("--force", is_flag=True, help="Re-provision even if envs are already present and current.")
def init_machine(processing_repo: str, device: str, manager: str | None,
                 container_runtime: str | None, singularity_module: str | None,
                 gpu_image: str | None, gpu_image_uri: str | None, conda_env: str | None,
                 gpu_conda_env: str | None, policy: str, install_processing: bool,
                 force: bool) -> None:
    """Make THIS machine ready for MuAgene in one operator-facing command.

    Bootstrap entry point: probes capabilities, writes ~/.muagene/machine.config,
    provisions the CPU env from the committed lock (NO science site.config needed),
    installs both agent packages into it, pulls the GPU image if requested, validates,
    and prints a readiness report. This is the FIRST thing run on a fresh machine —
    Execution-MuAgent owns non-scientific infrastructure end to end.

    Operator-facing by design: unlike the run-time lifecycle commands (which write
    findings into internal/hpc_monitor/ and never address the user), init-machine
    prints structured results to stdout — there is no Processing agent or run dir yet.
    """
    from . import environment, machine
    cfg = machine.detect_machine_config(
        processing_repo, manager=manager, container_runtime=container_runtime,
        singularity_module=singularity_module, gpu_image=gpu_image,
        gpu_image_uri=gpu_image_uri, policy=policy, conda_env=conda_env,
        gpu_conda_env=gpu_conda_env)
    mc_path = machine.write_machine_config(cfg)
    click.echo(f"machine profile: {mc_path}")
    click.echo(f"  manager={cfg.manager} runtime={cfg.container_runtime} "
               f"scheduler={cfg.scheduler} gpu_present={cfg.gpu_present}")
    if cfg.manager is None:
        raise click.ClickException(
            "No conda env manager (micromamba/mamba/conda) on PATH; cannot provision. "
            "Install one (e.g. miniforge) and re-run.")

    devices = ["cpu", "gpu"] if device == "both" else [device]
    if "gpu" in devices and not cfg.gpu_image_uri:
        raise click.ClickException(
            "device includes gpu but no --gpu-image-uri given. The GPU image is PULLED "
            "from a pinned, centrally-published registry reference — no machine builds it "
            "locally. Pass --gpu-image-uri docker://<registry>/muagene-gpu:<tag>.")

    failed: list[str] = []
    for dev in devices:
        site_cfg = machine.synthesize_site_config(processing_repo, cfg, device=dev)
        spec = environment.resolve_env_spec(site_cfg, processing_repo, dev)
        click.echo(f"[{dev}] provider={spec.provider} target={spec.image or spec.env_name}")
        res = environment.provision_env(spec, site_cfg, force=force)
        click.echo(f"[{dev}] {res.get('action', 'noop')} -> {res.get('status')}")
        if res.get("status") == "failed":
            click.echo((res.get("stderr") or "")[-1200:], err=True)
            failed.append(dev)
            continue
        val = environment.validate_env(spec, site_cfg)
        for f in val["findings"]:
            click.echo(f"[{dev}] {f['severity']}: {f['message']}", err=(f["severity"] == "error"))
        if not val["ok"]:
            failed.append(dev)

    # Layer the agent packages onto the CPU env so `Processing-MuAgent submit` (which
    # spawns `python -m execution_muagent`) runs both CLIs from the science env.
    if install_processing and "cpu" in devices and "cpu" not in failed:
        cpu_cfg = machine.synthesize_site_config(processing_repo, cfg, device="cpu")
        env_name = environment.resolve_env_spec(cpu_cfg, processing_repo, "cpu").env_name
        for label, repo in (("Processing-MuAgent", processing_repo),
                            ("Execution-MuAgent", machine.EXECUTION_REPO_ROOT)):
            pres = environment.pip_install_editable(cfg.manager, env_name, repo)
            click.echo(f"[cpu] pip install -e {label} -> rc={pres.get('returncode')}")
            if pres.get("returncode") != 0:
                click.echo((pres.get("stderr") or "")[-800:], err=True)
                failed.append(f"pip:{label}")

    if failed:
        raise click.ClickException(
            "init-machine did not complete cleanly for: " + ", ".join(failed))
    click.echo("Machine ready.")


if __name__ == "__main__":
    main()
