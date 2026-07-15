"""Runtime/status reporting — render `status` and `hpc-status` for the agent.

Both are read-only windows onto run state: `status` derives per-step pipeline state from
Snakemake logs + approval sentinels; `hpc-status` additionally reads Execution-MuAgent's
latest_snapshot.json (the machine contract) and prints the report-and-repoll fingerprint
+ single next-action line. Neither mutates run state.
"""
from __future__ import annotations

import os
import sys
import time

import click

from . import hpc, stage_progress as _sp
from .run_paths import RunPaths


def run_status(paths: RunPaths, *, watch: bool, interval: float) -> None:
    """Print per-step pipeline state once, or (with watch) poll until something changes."""
    def _print(states: list[tuple[str, str, str]]) -> None:
        for label, task, st in states:
            click.echo(f"  {label:18s}  {task:30s}  {st}")

    if not watch:
        _print(_sp.stage_states(paths))
        return

    sys.stdout.reconfigure(line_buffering=True)
    last: list[tuple[str, str, str]] | None = None
    while True:
        states = _sp.stage_states(paths)
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
            click.echo("\n→ a step was cancelled by the HPC monitor; run `hpc-status` "
                       "for the structured finding, confirmed-dead reason, and recovery "
                       "action. Re-`submit` after resolving the finding.")
            return
        if any(st == "awaiting_approval" for _, _, st in states):
            click.echo("\n→ a review gate is awaiting approval; review deliverables and run "
                       "`Processing-MuAgent approve <stage>` (e.g. qc_review).")
            return
        if paths.run_manifest_json.exists():
            click.echo("\n→ run_manifest.json present; pipeline complete.")
            return
        time.sleep(max(2.0, interval))


def run_hpc_status(paths: RunPaths) -> None:
    """Report HPC job health, monitor findings, and per-step pipeline state (one-shot).

    Reads only structured JSON (latest_snapshot.json + latest_submission.json) and prints
    once. Drives report-and-repoll: prints a stable `State:` fingerprint and the single
    next-action token the agent keys on. No poll loop — the daemon does the monitoring.
    """
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

    states = _sp.stage_states(paths)
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
