"""HPC job registration, monitoring, and diagnostics.

The monitor is intentionally scheduler-light: it uses the scheduler only for
state, and uses filesystem progress as the main signal that a workflow is alive.
All scheduler calls are bounded by short subprocess timeouts so the monitor
cannot hang behind a stuck `squeue`/`qstat` call.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


TERMINAL_ERROR_MARKERS = (
    "MissingInputException",
    "MissingOutputException",
    "LockException",
    "WorkflowError",
    "Error in rule",
    "Traceback",
    "OSError",
    "FileNotFoundError",
    "snakemake exited with",
    "DUE TO TIME LIMIT",
    "CANCELLED",
    "OUT_OF_MEMORY",
    "out of memory",
    "oom",
    "Killed",
)

SCHEDULER_FAILED_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "SPECIAL_EXIT",
    "TIMEOUT",
}

SCHEDULER_RUNNING_STATES = {
    "CONFIGURING",
    "COMPLETING",
    "PENDING",
    "RESIZING",
    "RUNNING",
    "SUSPENDED",
}


@dataclass
class Submission:
    agent: str
    executor: str
    job_id: str
    run_dir: str
    config: str
    target: str
    repo_root: str
    log_path: str
    submitted_at: str


@dataclass
class MonitorFinding:
    severity: str
    code: str
    message: str


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def run_monitor_dir(run_dir: Path | str) -> Path:
    return Path(run_dir) / "internal" / "hpc_monitor"


def register_submission(submission: Submission) -> Path:
    """Record an HPC submission in the run-local monitor registry."""
    monitor_dir = run_monitor_dir(submission.run_dir)
    monitor_dir.mkdir(parents=True, exist_ok=True)
    record_path = monitor_dir / "submissions.jsonl"
    with record_path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(asdict(submission), sort_keys=True) + "\n")
    latest_path = monitor_dir / "latest_submission.json"
    latest_path.write_text(json.dumps(asdict(submission), indent=2, sort_keys=True) + "\n")
    return record_path


def load_latest_submission(run_dir: Path | str, job_id: str | None = None) -> Submission:
    monitor_dir = run_monitor_dir(run_dir)
    records_path = monitor_dir / "submissions.jsonl"
    if not records_path.exists():
        raise FileNotFoundError(f"No HPC submissions registered under {monitor_dir}")
    selected: dict[str, Any] | None = None
    for line in records_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        record = json.loads(line)
        if job_id is None or str(record.get("job_id")) == str(job_id):
            selected = record
    if selected is None:
        raise FileNotFoundError(f"No registered submission for job {job_id!r} in {records_path}")
    return Submission(**selected)


def _run_cmd(cmd: list[str], timeout_s: int) -> tuple[int | None, str, str]:
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s, check=False)
        return result.returncode, result.stdout.strip(), result.stderr.strip()
    except subprocess.TimeoutExpired:
        return None, "", f"timeout after {timeout_s}s: {' '.join(cmd)}"
    except OSError as exc:
        return 127, "", str(exc)


def _normalize_job_id(job_id: str) -> str:
    return str(job_id).split(";", 1)[0].split(".", 1)[0].strip()


def query_slurm(job_id: str, timeout_s: int = 5) -> dict[str, str | None]:
    jid = _normalize_job_id(job_id)
    info: dict[str, str | None] = {
        "scheduler": "slurm",
        "job_id": jid,
        "state": None,
        "elapsed": None,
        "timelimit": None,
        "exit_code": None,
        "reason": None,
        "query_error": None,
    }
    if shutil.which("sacct"):
        rc, out, err = _run_cmd(
            ["sacct", "-j", jid, "-X", "-n", "-P", "-o", "State,Elapsed,Timelimit,ExitCode"],
            timeout_s,
        )
        if rc is None:
            info["query_error"] = err
        elif rc == 0 and out:
            parts = out.splitlines()[0].split("|")
            info["state"] = parts[0] if len(parts) > 0 else None
            info["elapsed"] = parts[1] if len(parts) > 1 else None
            info["timelimit"] = parts[2] if len(parts) > 2 else None
            info["exit_code"] = parts[3] if len(parts) > 3 else None
            return info
    if shutil.which("squeue"):
        rc, out, err = _run_cmd(["squeue", "-j", jid, "-h", "-o", "%T|%M|%l|%R"], timeout_s)
        if rc is None:
            info["query_error"] = err
        elif rc == 0 and out:
            parts = out.splitlines()[0].split("|")
            info["state"] = parts[0] if len(parts) > 0 else None
            info["elapsed"] = parts[1] if len(parts) > 1 else None
            info["timelimit"] = parts[2] if len(parts) > 2 else None
            info["reason"] = parts[3] if len(parts) > 3 else None
    return info


def query_pbs(job_id: str, timeout_s: int = 5) -> dict[str, str | None]:
    info: dict[str, str | None] = {
        "scheduler": "pbs",
        "job_id": job_id,
        "state": None,
        "elapsed": None,
        "timelimit": None,
        "exit_code": None,
        "reason": None,
        "query_error": None,
    }
    if not shutil.which("qstat"):
        info["query_error"] = "qstat not found"
        return info
    rc, out, err = _run_cmd(["qstat", "-f", job_id], timeout_s)
    if rc is None:
        info["query_error"] = err
        return info
    if rc != 0:
        rc, out, err = _run_cmd(["qstat", "-fx", job_id], timeout_s)
        if rc is None:
            info["query_error"] = err
            return info
    if out:
        for line in out.splitlines():
            if "job_state" in line and "=" in line:
                info["state"] = line.split("=", 1)[1].strip()
            elif "Exit_status" in line and "=" in line:
                info["exit_code"] = line.split("=", 1)[1].strip()
            elif "resources_used.walltime" in line and "=" in line:
                info["elapsed"] = line.split("=", 1)[1].strip()
            elif "Resource_List.walltime" in line and "=" in line:
                info["timelimit"] = line.split("=", 1)[1].strip()
    return info


def query_scheduler(executor: str, job_id: str, timeout_s: int = 5) -> dict[str, str | None]:
    if executor == "slurm":
        return query_slurm(job_id, timeout_s)
    if executor == "pbs":
        return query_pbs(job_id, timeout_s)
    return {"scheduler": executor, "job_id": job_id, "state": None, "query_error": "unsupported executor"}


def _path_status(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"path": str(path), "exists": False}
    stat = path.stat()
    return {
        "path": str(path),
        "exists": True,
        "size": stat.st_size,
        "mtime": stat.st_mtime,
        "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def _newest_file(paths: list[Path]) -> Path | None:
    existing = [p for p in paths if p.exists() and p.is_file()]
    if not existing:
        return None
    return max(existing, key=lambda p: p.stat().st_mtime)


def _submitted_epoch(submission: Submission) -> float:
    try:
        return datetime.strptime(submission.submitted_at, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc).timestamp()
    except ValueError:
        return 0.0


def snakemake_log_roots(submission: Submission) -> list[Path]:
    """Return run-local log roots first, with legacy repo-level roots as fallback."""
    run_dir = Path(submission.run_dir)
    repo_root = Path(submission.repo_root)
    return [
        run_dir / "internal" / "snakemake" / ".snakemake" / "log",
        run_dir / "internal" / "snakemake" / ".snakemake" / "slurm_logs",
        repo_root / ".snakemake" / "log",
        repo_root / ".snakemake" / "slurm_logs",
        repo_root / "logs",
    ]


def _file_mentions(path: Path, needle: str, max_bytes: int = 2_000_000) -> bool:
    try:
        with path.open("rb") as fh:
            data = fh.read(max_bytes)
        return needle in data.decode("utf-8", errors="replace")
    except OSError:
        return False


def discover_log_files(submission: Submission) -> list[Path]:
    submitted_after = max(0.0, _submitted_epoch(submission) - 60.0)
    run_dir_text = str(Path(submission.run_dir))
    files: list[Path] = []
    head_log = Path(submission.log_path)
    if head_log.exists() and head_log.is_file():
        files.append(head_log)
    roots = snakemake_log_roots(submission)
    run_local_roots = {
        Path(submission.run_dir) / "internal" / "snakemake" / ".snakemake" / "log",
        Path(submission.run_dir) / "internal" / "snakemake" / ".snakemake" / "slurm_logs",
    }
    for root in roots:
        if root.exists():
            for p in root.rglob("*"):
                if not p.is_file() or p.stat().st_mtime < submitted_after:
                    continue
                if root in run_local_roots or _file_mentions(p, run_dir_text):
                    files.append(p)
    return sorted(set(files), key=lambda p: str(p))


def discover_progress_files(submission: Submission) -> list[Path]:
    run_dir = Path(submission.run_dir)
    files = [
        Path(submission.log_path),
        run_dir / "internal" / "log.jsonl",
        run_dir / "internal" / "parameters.yaml",
    ]
    files.extend(discover_log_files(submission))
    for root in (run_dir / "internal" / "artifacts", run_dir / "deliverables"):
        if root.exists():
            files.extend(p for p in root.rglob("*") if p.is_file())
    return files


def _read_tail(path: Path, max_bytes: int = 8000) -> str:
    if not path.exists() or not path.is_file():
        return ""
    with path.open("rb") as fh:
        if path.stat().st_size > max_bytes:
            fh.seek(-max_bytes, os.SEEK_END)
        data = fh.read()
    return data.decode("utf-8", errors="replace")


def _extract_markers(text: str) -> list[str]:
    found: list[str] = []
    for marker in TERMINAL_ERROR_MARKERS:
        if marker.lower() in text.lower():
            found.append(marker)
    return found


def discover_child_job_ids(submission: Submission) -> list[str]:
    """Extract child scheduler ids from all known head/Snakemake logs."""
    head_id = _normalize_job_id(submission.job_id)
    ids: set[str] = set()
    for path in discover_log_files(submission):
        for job_id in parse_job_ids_from_log(path):
            norm = _normalize_job_id(job_id)
            if norm and norm != head_id:
                ids.add(norm)
    return sorted(ids)


def collect_snapshot(submission: Submission, scheduler_timeout_s: int = 5) -> dict[str, Any]:
    scheduler = query_scheduler(submission.executor, submission.job_id, scheduler_timeout_s)
    progress_files = discover_progress_files(submission)
    newest = _newest_file(progress_files)
    head_log = Path(submission.log_path)
    head_tail = _read_tail(head_log)
    latest_snakemake_log = _newest_file(
        [p for p in discover_log_files(submission) if ".snakemake" in p.parts and p.suffix == ".log"]
    )
    latest_snakemake_tail = _read_tail(latest_snakemake_log) if latest_snakemake_log else ""
    markers = sorted(set(_extract_markers(head_tail) + _extract_markers(latest_snakemake_tail)))
    child_job_ids = discover_child_job_ids(submission)
    writing_files = [
        _path_status(p)
        for p in (Path(submission.run_dir) / "internal" / "artifacts").rglob("*.writing")
    ] if (Path(submission.run_dir) / "internal" / "artifacts").exists() else []
    return {
        "submission": asdict(submission),
        "checked_at": utc_now(),
        "scheduler": scheduler,
        "head_log": _path_status(head_log),
        "latest_progress_file": _path_status(newest) if newest else None,
        "latest_snakemake_log": _path_status(latest_snakemake_log) if latest_snakemake_log else None,
        "child_job_ids": child_job_ids,
        "error_markers": markers,
        "writing_files": writing_files,
        "head_log_tail": head_tail,
        "latest_snakemake_tail": latest_snakemake_tail,
    }


def evaluate_snapshot(snapshot: dict[str, Any], stale_minutes: float) -> list[MonitorFinding]:
    findings: list[MonitorFinding] = []
    scheduler = snapshot.get("scheduler") or {}
    state_parts = str(scheduler.get("state") or "").split()
    state = state_parts[0] if state_parts else ""
    query_error = scheduler.get("query_error")
    if query_error:
        findings.append(MonitorFinding("warning", "scheduler_query_failed", str(query_error)))
    if state in SCHEDULER_FAILED_STATES:
        findings.append(MonitorFinding("error", "scheduler_failed", f"Scheduler state is {state}."))
    if snapshot.get("error_markers"):
        findings.append(
            MonitorFinding(
                "error",
                "workflow_error_marker",
                "Workflow log contains: " + ", ".join(snapshot["error_markers"]),
            )
        )
    latest = snapshot.get("latest_progress_file") or {}
    mtime = latest.get("mtime")
    if mtime:
        age_min = (time.time() - float(mtime)) / 60.0
        if state in SCHEDULER_RUNNING_STATES and age_min >= stale_minutes:
            findings.append(
                MonitorFinding(
                    "error",
                    "no_filesystem_progress",
                    f"No monitored file changed for {age_min:.1f} min while scheduler state is {state}.",
                )
            )
        elif not state and age_min >= stale_minutes:
            findings.append(
                MonitorFinding(
                    "error",
                    "unknown_scheduler_state_no_progress",
                    f"Scheduler state is unknown and no monitored file changed for {age_min:.1f} min.",
                )
            )
    elif state in SCHEDULER_RUNNING_STATES:
        findings.append(MonitorFinding("warning", "no_progress_files", "No progress files found yet."))
    elif not state:
        age_min = (time.time() - _submitted_epoch(Submission(**snapshot["submission"]))) / 60.0
        if age_min >= stale_minutes:
            findings.append(
                MonitorFinding(
                    "error",
                    "unknown_scheduler_state_no_progress_files",
                    f"Scheduler state is unknown and no progress files appeared for {age_min:.1f} min.",
                )
            )
    if state == "COMPLETING":
        findings.append(
            MonitorFinding(
                "warning",
                "scheduler_completing",
                "Job is in COMPLETING; on this cluster that can indicate epilog or NFS cleanup stall.",
            )
        )
    return findings


def cancel_job(executor: str, job_id: str, timeout_s: int = 5) -> dict[str, Any]:
    if executor == "slurm":
        cmd = ["scancel", _normalize_job_id(job_id)]
    elif executor == "pbs":
        cmd = ["qdel", job_id]
    else:
        return {"attempted": False, "reason": f"unsupported executor {executor!r}"}
    rc, out, err = _run_cmd(cmd, timeout_s)
    return {"attempted": True, "returncode": rc, "stdout": out, "stderr": err, "cmd": cmd}


def cancel_submission_jobs(submission: Submission, child_job_ids: list[str], timeout_s: int = 5) -> dict[str, Any]:
    """Cancel child jobs first, then the registered head job."""
    actions: list[dict[str, Any]] = []
    for job_id in child_job_ids:
        actions.append({"job_id": job_id, "role": "child", **cancel_job(submission.executor, job_id, timeout_s)})
    actions.append({
        "job_id": submission.job_id,
        "role": "head",
        **cancel_job(submission.executor, submission.job_id, timeout_s),
    })
    return {
        "attempted": True,
        "executor": submission.executor,
        "head_job_id": submission.job_id,
        "child_job_ids": child_job_ids,
        "actions": actions,
    }


def render_report(snapshot: dict[str, Any], findings: list[MonitorFinding], cancel_result: dict[str, Any] | None) -> str:
    sub = snapshot["submission"]
    sched = snapshot["scheduler"]
    lines = [
        "# HPC Monitor Report",
        "",
        f"- Run: `{sub['run_dir']}`",
        f"- Agent: `{sub['agent']}`",
        f"- Executor: `{sub['executor']}`",
        f"- Head job: `{sub['job_id']}`",
        f"- Target: `{sub['target']}`",
        f"- Checked: `{snapshot['checked_at']}`",
        f"- Scheduler state: `{sched.get('state') or 'unknown'}`",
        f"- Elapsed / limit: `{sched.get('elapsed') or 'unknown'}` / `{sched.get('timelimit') or 'unknown'}`",
        "",
        "## Findings",
        "",
    ]
    if findings:
        for finding in findings:
            lines.append(f"- **{finding.severity.upper()} `{finding.code}`**: {finding.message}")
    else:
        lines.append("- No monitor problems detected in this check.")
    if cancel_result is not None:
        lines.extend(["", "## Kill Action", "", f"```json\n{json.dumps(cancel_result, indent=2)}\n```"])
    lines.extend(
        [
            "",
            "## Progress",
            "",
            f"- Head log: `{snapshot['head_log'].get('path')}`",
            f"- Latest progress file: `{(snapshot.get('latest_progress_file') or {}).get('path', 'none')}`",
            f"- Latest Snakemake log: `{(snapshot.get('latest_snakemake_log') or {}).get('path', 'none')}`",
            f"- Child job ids: `{', '.join(snapshot.get('child_job_ids') or []) or 'none'}`",
            "",
            "## Last Head Log Lines",
            "",
            "```text",
            _tail_lines(snapshot.get("head_log_tail") or "", 80),
            "```",
        ]
    )
    if snapshot.get("latest_snakemake_tail"):
        lines.extend(
            [
                "",
                "## Last Snakemake Log Lines",
                "",
                "```text",
                _tail_lines(snapshot["latest_snakemake_tail"], 80),
                "```",
            ]
        )
    return "\n".join(lines) + "\n"


def _tail_lines(text: str, n: int) -> str:
    lines = text.splitlines()
    return "\n".join(lines[-n:])


def write_report(
    submission: Submission,
    snapshot: dict[str, Any],
    findings: list[MonitorFinding],
    cancel_result: dict[str, Any] | None = None,
) -> Path:
    report_dir = run_monitor_dir(submission.run_dir) / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    report_path = report_dir / f"{submission.job_id}_{stamp}.md"
    text = render_report(snapshot, findings, cancel_result)
    report_path.write_text(text, encoding="utf-8")
    latest = run_monitor_dir(submission.run_dir) / "latest_report.md"
    latest.write_text(text, encoding="utf-8")
    state_path = run_monitor_dir(submission.run_dir) / "latest_snapshot.json"
    state_path.write_text(json.dumps(snapshot, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report_path


def monitor_once(
    submission: Submission,
    *,
    stale_minutes: float,
    scheduler_timeout_s: int,
    kill_on_hang: bool,
) -> tuple[dict[str, Any], list[MonitorFinding], dict[str, Any] | None, Path | None]:
    snapshot = collect_snapshot(submission, scheduler_timeout_s)
    findings = evaluate_snapshot(snapshot, stale_minutes)
    killable = any(
        f.code in {
            "scheduler_failed",
            "no_filesystem_progress",
            "unknown_scheduler_state_no_progress",
            "unknown_scheduler_state_no_progress_files",
        }
        for f in findings
    )
    child_job_ids = snapshot.get("child_job_ids") or []
    cancel_result = (
        cancel_submission_jobs(submission, child_job_ids, scheduler_timeout_s)
        if kill_on_hang and killable else None
    )
    report_path = write_report(submission, snapshot, findings, cancel_result) if findings or cancel_result else None
    return snapshot, findings, cancel_result, report_path


def monitor_watch(
    submission: Submission,
    *,
    interval_s: float,
    stale_minutes: float,
    scheduler_timeout_s: int,
    kill_on_hang: bool,
    max_checks: int | None,
) -> Path | None:
    checks = 0
    last_report: Path | None = None
    while True:
        snapshot, findings, cancel_result, report_path = monitor_once(
            submission,
            stale_minutes=stale_minutes,
            scheduler_timeout_s=scheduler_timeout_s,
            kill_on_hang=kill_on_hang,
        )
        if report_path is not None:
            last_report = report_path
        if cancel_result is not None or any(f.severity == "error" for f in findings):
            return last_report
        state_parts = str((snapshot.get("scheduler") or {}).get("state") or "").split()
        state = state_parts[0] if state_parts else ""
        if state == "COMPLETED":
            return last_report
        checks += 1
        if max_checks is not None and checks >= max_checks:
            return last_report
        time.sleep(max(5.0, interval_s))


def parse_job_ids_from_log(path: Path | str) -> list[str]:
    """Extract scheduler ids from a Snakemake/head-job log."""
    p = Path(path)
    if not p.exists():
        return []
    text = p.read_text(encoding="utf-8", errors="replace")
    patterns = [
        r"SLURM jobid ([0-9]+)",
        r"Submitted batch job ([0-9]+)",
        r"Submitted .*head-job: ([0-9A-Za-z_.-]+)",
        r"Submitted .*: ([0-9A-Za-z_.-]+)",
    ]
    ids: list[str] = []
    for pat in patterns:
        ids.extend(re.findall(pat, text))
    return sorted(set(ids))
