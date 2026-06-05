"""HPC helpers — profile paths, head-job submission, scheduler detection.

Used by `executor.cli` for the `--executor` flag and the `submit` command.
Keeps cluster knowledge out of the CLI module itself.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Literal


Executor = Literal["local", "pbs", "slurm"]

# Repo root: Processing-MuAgent/ — derived from this file's location.
REPO_ROOT: Path = Path(__file__).resolve().parent.parent

PROFILE_DIR = {
    "pbs":   REPO_ROOT / "workflow" / "profiles" / "pbs",
    "slurm": REPO_ROOT / "workflow" / "profiles" / "slurm",
}

RUNNER_SCRIPT = {
    "pbs":   REPO_ROOT / "scripts" / "runner.pbs",
    "slurm": REPO_ROOT / "scripts" / "runner.slurm",
}

LAUNCHER = REPO_ROOT / "scripts" / "launch_runner.sh"

# Default directory for head-job and PBS child-job stdout/stderr (see pbs-submit.sh).
DEFAULT_LOG_DIR = REPO_ROOT / "logs"

# Snakemake 9 defaults --shared-fs-usage to ALL, including storage-local-copies.
# On shared NFS (runs/... under /home/.../mnt/storage/...) that makes cluster
# child jobs run in --mode remote and enter a post-job "Storing output in
# storage." phase. Our stages already write declared outputs directly to NFS
# (see executor/io.write_h5ad_safe), so the storage sync can hang indefinitely
# while SLURM still shows the child as RUNNING. Omit storage-local-copies and
# keep input-output on the shared filesystem. See workflow/profiles/*/config.yaml.
SNAKEMAKE_SHARED_FS_USAGE: tuple[str, ...] = (
    "persistence",
    "input-output",
    "software-deployment",
    "software-deployment-cache",
    "sources",
    "source-cache",
)


def snakemake_cluster_cli_args() -> list[str]:
    """Extra snakemake CLI flags for PBS/SLURM orchestration on shared NFS."""
    return ["--shared-fs-usage", *SNAKEMAKE_SHARED_FS_USAGE]


def sanitize_snakemake_jobscript_text(text: str) -> str:
    """Rewrite Snakemake job scripts for safe, efficient execution on shared NFS clusters.

    Removes `storage-local-copies` and `--local-storage-prefix` options, and replaces
    `--mode remote` with `--mode subprocess`. This avoids unnecessary local storage
    use and prevents post-job hangs during output syncing, as all I/O should occur
    directly over NFS.
    """
    text = re.sub(r"(?<=\s)storage-local-copies(?=\s)", "", text)
    text = re.sub(r"--mode\s+'remote'", "--mode 'subprocess'", text)
    text = re.sub(r"--mode\s+remote(?=\s)", "--mode subprocess ", text)
    text = re.sub(r"\s--local-storage-prefix\s+\S+", "", text)
    return text


def sanitize_snakemake_jobscript(path: Path | str) -> bool:
    """Sanitize a generated Snakemake child jobscript in-place.

    Returns True when the file was changed.
    """
    p = Path(path)
    text = p.read_text()
    sanitized = sanitize_snakemake_jobscript_text(text)
    if sanitized == text:
        return False
    p.write_text(sanitized)
    return True


def snakemake_lock_files(workdir: Path | str) -> list[Path]:
    """Return existing Snakemake lock files for a run workdir."""
    lock_dir = Path(workdir) / ".snakemake" / "locks"
    if not lock_dir.exists():
        return []
    return sorted(p for p in lock_dir.iterdir() if p.is_file())


def snakemake_processes_for_workdir(workdir: Path | str) -> list[tuple[int, str]]:
    """Best-effort local process check for Snakemake using a run workdir."""
    workdir_s = str(Path(workdir).resolve())
    try:
        result = subprocess.run(
            ["ps", "-eo", "pid=,args="],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return []
    if result.returncode != 0:
        return []
    current = os.getpid()
    matches: list[tuple[int, str]] = []
    for line in result.stdout.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        pid_s, _, args = stripped.partition(" ")
        try:
            pid = int(pid_s)
        except ValueError:
            continue
        if pid == current:
            continue
        is_snakemake_cmd = (
            " -m snakemake" in args
            or args.startswith("snakemake ")
            or "/snakemake " in args
        )
        if is_snakemake_cmd and workdir_s in args:
            matches.append((pid, args))
    return matches


def resolve_log_dir() -> Path:
    """Return the scheduler log directory, creating it if needed."""
    raw = os.environ.get("PMA_LOG_DIR", "logs")
    log_dir = Path(raw)
    if not log_dir.is_absolute():
        log_dir = REPO_ROOT / log_dir
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir


def head_job_log_path(executor: Executor) -> Path:
    """Default head-job log path for the given scheduler."""
    log_dir = resolve_log_dir()
    if executor == "slurm":
        return log_dir / "pma_runner-%j.out"
    return log_dir


def profile_path(executor: Executor) -> Path:
    """Return the snakemake profile directory for the given executor."""
    if executor == "local":
        raise ValueError("profile_path is not applicable for local executor")
    p = PROFILE_DIR[executor]
    if not p.exists():
        raise FileNotFoundError(
            f"snakemake profile dir not found: {p}. Run from a clean checkout.")
    return p


def detect_scheduler() -> Executor:
    """Best-effort detection of which scheduler is available on PATH.

    Returns 'pbs' if qsub is present, 'slurm' if sbatch is present, 'local' otherwise.
    Used for friendlier default behaviour when --executor is omitted on a known cluster.
    """
    if shutil.which("qsub"):
        return "pbs"
    if shutil.which("sbatch"):
        return "slurm"
    return "local"


def submitted_log_path(executor: Executor, output_log: Path | str, job_id: str) -> Path:
    """Return the concrete scheduler log path after scheduler id substitution."""
    path = Path(output_log)
    if not path.is_absolute():
        path = REPO_ROOT / path
    if executor == "slurm":
        rendered = str(path).replace("%j", str(job_id)).replace("%A", str(job_id))
        return Path(rendered)
    return path


def kill_existing_supervisor(run_dir: Path | str) -> bool:
    """Kill the running supervisor daemon for a run, if one exists.

    Reads internal/hpc_monitor/monitor.pid; sends SIGTERM if the process is alive
    and removes the PID file. Returns True when a live supervisor was killed.
    Gate-resubmit cycles are safe: the old supervisor exits when the job finishes,
    so os.kill(pid, 0) raises OSError (dead PID) and nothing is killed.
    """
    import signal as _signal
    pid_path = Path(run_dir) / "internal" / "hpc_monitor" / "monitor.pid"
    if not pid_path.exists():
        return False
    try:
        pid = int(pid_path.read_text().strip())
        os.kill(pid, 0)           # raises OSError if dead
        os.kill(pid, _signal.SIGTERM)
        pid_path.unlink(missing_ok=True)
        return True
    except (ValueError, OSError):
        return False


def _wait_for_manifest(run_dir: Path | str, timeout_s: float = 90.0) -> str | None:
    """Poll execution_manifest.jsonl until a job_id appears or timeout_s elapses.

    Returns the job_id string, or None on timeout. The caller must warn explicitly
    on None — do NOT silently log "unknown".
    """
    import time as _time
    deadline = _time.monotonic() + timeout_s
    while _time.monotonic() < deadline:
        job_id = last_manifest_job_id(run_dir)
        if job_id:
            return job_id
        _time.sleep(0.5)
    return None


def start_supervisor_daemon(
    run_dir: Path | str,
    cmd: list[str],
    env: dict[str, str],
) -> dict[str, object] | None:
    """Launch execute-spec (or resume-monitor) as a detached background daemon.

    Uses a timestamped log file with a `monitor.log` symlink pointing to the latest,
    so history is preserved across gate-resubmit cycles. `start_new_session=True`
    puts the child in a new OS session with no controlling terminal:
    - SSH SIGHUP only reaches the old session → daemon survives SSH disconnect.
    - Ctrl-C SIGINT targets the foreground process group only → daemon unaffected.
    Caveat: sites with KillUserProcesses=yes (systemd) kill all user processes on
    logout regardless of session; use tmux/screen on those systems.

    Returns {"pid": int, "log": str} on success, None on failure.
    """
    from datetime import datetime as _dt
    log_dir = Path(run_dir) / "internal" / "hpc_monitor"
    log_dir.mkdir(parents=True, exist_ok=True)
    ts = _dt.now().strftime("%Y%m%d_%H%M%S")
    log_path = log_dir / f"monitor_{ts}.log"
    symlink = log_dir / "monitor.log"
    if symlink.exists() or symlink.is_symlink():
        symlink.unlink()
    symlink.symlink_to(log_path.name)
    pid_path = log_dir / "monitor.pid"
    try:
        with open(log_path, "w") as log_fh:
            proc = subprocess.Popen(
                cmd, env=env,
                stdout=log_fh, stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL,
                start_new_session=True,
            )
        pid_path.write_text(str(proc.pid))
        return {"pid": proc.pid, "log": str(log_path)}
    except (OSError, subprocess.SubprocessError):
        return None


def _execution_muagent_env() -> dict[str, str] | None:
    """Return an environment that can import sibling Execution-MuAgent, if present."""
    exec_root = REPO_ROOT.parent / "Execution-MuAgent"
    if not (exec_root / "execution_muagent").exists():
        return None
    env = os.environ.copy()
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = str(exec_root) if not existing else f"{exec_root}:{existing}"
    return env


def last_manifest_job_id(run_dir: Path | str) -> str | None:
    """Return the most-recently recorded job_id from execution_manifest.jsonl, or None."""
    import json as _json
    p = Path(run_dir) / "internal" / "hpc_monitor" / "execution_manifest.jsonl"
    if not p.exists():
        return None
    for line in reversed(p.read_text().splitlines()):
        line = line.strip()
        if not line:
            continue
        try:
            return _json.loads(line).get("job_id")
        except Exception:
            continue
    return None


def submit_via_execution_muagent(
    spec_path: Path | str,
    site_config_path: Path | str,
    run_dir: Path | str,
    target: str,
    *,
    watch: bool = False,
    kill_on_hang: bool = True,
) -> dict[str, object] | None:
    """Delegate HPC submission to Execution-MuAgent via execute-spec.

    watch=False: captures stdout/stderr for job-ID parsing; returns within 120 s.
    No supervisor daemon is started — no stall detection, no auto-cancel.

    watch=True: starts execute-spec --watch as a detached daemon (start_new_session=True).
    The daemon survives SSH disconnect (unless the site uses KillUserProcesses=yes;
    use tmux/screen there). Polls execution_manifest.jsonl for up to 90 s to read the
    job_id, then returns. Execution-MuAgent output goes to monitor_<ts>.log — it is
    NOT streamed to the terminal (Processing-MuAgent is the sole UX layer).

    Returns None when Execution-MuAgent is unavailable or on error.
    The caller must raise a hard error on None — there is no direct-submit fallback.
    """
    env = _execution_muagent_env()
    if env is None:
        return None
    cmd = [
        sys.executable, "-m", "execution_muagent.cli", "execute-spec",
        "--spec", str(Path(spec_path).resolve()),
        "--site-config", str(Path(site_config_path).resolve()),
        "--run-dir", str(Path(run_dir).resolve()),
        "--repo-root", str(REPO_ROOT),
        "--target", target,
    ]
    if watch:
        cmd.append("--watch")
    if kill_on_hang:
        cmd.append("--kill-on-hang")
    else:
        cmd.append("--no-kill-on-hang")
    try:
        if watch:
            killed = kill_existing_supervisor(run_dir)
            if killed:
                sys.stderr.write("[submit] Stopped previous supervisor before starting new one.\n")
            result = start_supervisor_daemon(run_dir, cmd, env)
            if result is None:
                return None
            job_id = _wait_for_manifest(run_dir, timeout_s=90.0)
            return {
                "watch": True,
                "pid": result["pid"],
                "log": result["log"],
                "job_id": job_id,   # None → caller warns explicitly
                "stdout": "",
                "stderr": "",
            }
        else:
            r = subprocess.run(cmd, env=env, capture_output=True, text=True,
                               timeout=120, check=True)
            return {"watch": False, "stdout": r.stdout, "stderr": r.stderr}
    except (subprocess.SubprocessError, OSError):
        return None


def env_diagnostics() -> dict[str, str | None]:
    """Snapshot of HPC-relevant env vars — used by `Processing-MuAgent hpc-info`
    to show the user what's wired up.
    """
    keys = (
        "PMA_PBS_QUEUE", "PMA_PBS_PROJECT",
        "PMA_SLURM_PARTITION", "PMA_SLURM_ACCOUNT",
        "PMA_RESOURCES_SCALE", "PMA_CONDA_ENV", "PMA_LOG_DIR",
    )
    return {k: os.environ.get(k) for k in keys}


def _run_cmd(cmd: list[str], *, timeout: int = 15) -> str:
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, check=False,
        )
        if result.returncode == 0:
            return result.stdout
    except (subprocess.TimeoutExpired, OSError):
        pass
    return ""


def _parse_pbs_queues(text: str) -> list[str]:
    queues: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("-") or line.startswith("Queue"):
            continue
        parts = line.split()
        if parts and parts[0] not in {"Queue", "---"}:
            queues.append(parts[0])
    return sorted(set(queues))


def _parse_slurm_partitions(text: str) -> list[str]:
    parts: list[str] = []
    for line in text.splitlines():
        name = line.strip().split()[0] if line.strip() else ""
        if name:
            parts.append(name.rstrip("*"))
    return sorted(set(parts))


def _parse_slurm_accounts(text: str) -> list[str]:
    accounts: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("Account"):
            continue
        acct = line.split("|")[0].split()[0] if "|" in line else line.split()[0]
        if acct and acct not in {"Account", "----------"}:
            accounts.append(acct)
    return sorted(set(accounts))


def discover_site() -> dict[str, object]:
    """Probe the login node for scheduler type, queues/partitions, and accounts.

    Best-effort — individual probes may return empty lists when the scheduler
    CLI is unavailable or the site restricts listing.
    """
    detected = detect_scheduler()
    user = os.environ.get("USER") or os.environ.get("LOGNAME") or ""
    info: dict[str, object] = {
        "detected_scheduler": detected,
        "user": user,
        "current_env": env_diagnostics(),
        "pbs": {"queues": [], "projects": [], "suggested_queue": None, "suggested_project": None},
        "slurm": {"partitions": [], "accounts": [], "suggested_partition": None, "suggested_account": None},
    }

    if detected == "pbs":
        pbs = info["pbs"]
        assert isinstance(pbs, dict)
        queue_text = _run_cmd(["qstat", "-Q"])
        queues = _parse_pbs_queues(queue_text)
        pbs["queues"] = queues

        # Site-specific hints already in the environment.
        for key in ("PMA_PBS_QUEUE", "PBS_DEFAULT_QUEUE", "PBS_QUEUE"):
            val = os.environ.get(key)
            if val:
                pbs["suggested_queue"] = val
                break
        if not pbs["suggested_queue"] and queues:
            pbs["suggested_queue"] = queues[0]

        for key in ("PMA_PBS_PROJECT", "PBS_PROJECT", "PBS_ACCOUNT"):
            val = os.environ.get(key)
            if val:
                pbs["projects"] = [val]
                pbs["suggested_project"] = val
                break

        # Recent jobs may expose a project code (-P).
        if user:
            recent = _run_cmd(["qstat", "-u", user, "-f"])
            for line in recent.splitlines():
                if "Project" in line:
                    proj = line.split("=", 1)[-1].strip()
                    if proj and proj not in pbs["projects"]:
                        pbs["projects"].append(proj)
            if not pbs["suggested_project"] and pbs["projects"]:
                pbs["suggested_project"] = pbs["projects"][0]

    elif detected == "slurm":
        slurm = info["slurm"]
        assert isinstance(slurm, dict)
        part_text = _run_cmd(["sinfo", "-h", "-o", "%P"])
        partitions = _parse_slurm_partitions(part_text)
        slurm["partitions"] = partitions

        for key in ("PMA_SLURM_PARTITION", "SLURM_PARTITION"):
            val = os.environ.get(key)
            if val:
                slurm["suggested_partition"] = val
                break
        if not slurm["suggested_partition"] and partitions:
            slurm["suggested_partition"] = partitions[0]

        if user:
            acct_text = _run_cmd([
                "sacctmgr", "show", "assoc", f"where=user={user}",
                "format=Account,Partition", "-P", "-n",
            ])
            accounts = _parse_slurm_accounts(acct_text)
            slurm["accounts"] = accounts

        for key in ("PMA_SLURM_ACCOUNT", "SLURM_ACCOUNT"):
            val = os.environ.get(key)
            if val:
                slurm["suggested_account"] = val
                break
        if not slurm["suggested_account"] and slurm["accounts"]:
            slurm["suggested_account"] = slurm["accounts"][0]

    return info


def write_site_config(path: Path | str, *, mode: Executor, settings: dict[str, str | None]) -> Path:
    """Write site.config — the YAML platform description consumed by Execution-MuAgent.

    Processing-MuAgent writes this from confirmed user input; Execution-MuAgent
    reads it to render submission scripts without scheduler knowledge baked in.
    """
    import yaml  # type: ignore[import]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    scheduler_section: dict[str, dict[str, str | None]] = {}
    if mode == "slurm":
        scheduler_section["slurm"] = {
            "partition": settings.get("slurm_partition"),
            "account": settings.get("slurm_account"),
            "qos": None,
        }
    elif mode == "pbs":
        scheduler_section["pbs"] = {
            "queue": settings.get("pbs_queue"),
            "project": settings.get("pbs_project"),
        }
    try:
        scale = int(float(settings["resources_scale"])) if settings.get("resources_scale") else 1
    except (ValueError, TypeError):
        scale = 1
    config: dict[str, object] = {
        "schema_version": "1",
        "scheduler": mode,
        **scheduler_section,
        "common": {
            "resources_scale": scale,
            "conda_env": settings.get("conda_env"),
            "container": None,
            "scratch": None,
        },
    }
    path.write_text(yaml.safe_dump(config, default_flow_style=False, sort_keys=False))
    return path


def load_site_config(path: Path | str) -> dict[str, object]:
    """Load site.config YAML; returns empty dict if the file does not exist."""
    import yaml  # type: ignore[import]
    p = Path(path)
    if not p.exists():
        return {}
    with p.open() as f:
        return yaml.safe_load(f) or {}


def write_hpc_env(path: Path | str, site_config_path: Path | str) -> Path:
    """Write a source-able shell snippet with PMA_* exports derived from site.config.

    Derives all values from site.config so the two files cannot drift — hpc.env
    is always a shell-variable projection of site.config, not a parallel source.
    """
    cfg = load_site_config(site_config_path)
    mode = cfg.get("scheduler", "local")
    common = cfg.get("common", {}) or {}
    slurm = cfg.get("slurm", {}) or {}
    pbs = cfg.get("pbs", {}) or {}

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Processing-MuAgent HPC settings — source before submit/run on cluster:",
        "#   source deliverables/pre_run/config/hpc.env",
        f"# execution mode: {mode}",
        "",
    ]
    exports: dict[str, str | None] = {
        "PMA_PBS_QUEUE":        pbs.get("queue"),
        "PMA_PBS_PROJECT":      pbs.get("project"),
        "PMA_SLURM_PARTITION":  slurm.get("partition"),
        "PMA_SLURM_ACCOUNT":    slurm.get("account"),
        "PMA_RESOURCES_SCALE":  str(common["resources_scale"]) if common.get("resources_scale") else None,
        "PMA_CONDA_ENV":        common.get("conda_env"),
    }
    for env_key, val in exports.items():
        if val:
            lines.append(f"export {env_key}={val!r}")
    path.write_text("\n".join(lines) + "\n")
    return path


def load_execution_mode(parameters_path: Path | str) -> Executor:
    """Read execution.mode from parameters.yaml; default local."""
    import yaml
    p = Path(parameters_path)
    if not p.exists():
        return "local"
    with p.open() as f:
        params = yaml.safe_load(f) or {}
    entry = params.get("execution.mode") or {}
    mode = entry.get("value") if isinstance(entry, dict) else entry
    if mode in ("local", "pbs", "slurm"):
        return mode
    return "local"


def parse_hpc_env(path: Path | str) -> dict[str, str]:
    """Parse export PMA_*=... lines from an hpc.env shell snippet."""
    p = Path(path)
    if not p.exists():
        return {}
    out: dict[str, str] = {}
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line.startswith("export PMA_"):
            continue
        body = line[len("export ") :]
        if "=" not in body:
            continue
        key, _, raw = body.partition("=")
        val = raw.strip().strip("'\"")
        out[key] = val
    return out


def load_execution_settings(run_dir: Path | str) -> dict[str, object]:
    """Execution mode + HPC settings recorded for this run (for plan review)."""
    from .run_paths import RunPaths

    paths = RunPaths(Path(run_dir))
    mode = load_execution_mode(paths.parameters_yaml)
    hpc_env_path = paths.hpc_env_sh
    from_file = parse_hpc_env(hpc_env_path) if hpc_env_path.exists() else {}
    live = env_diagnostics()

    def _get(env_key: str, field: str) -> str | None:
        return from_file.get(env_key) or live.get(env_key)

    settings: dict[str, str | None] = {
        "pbs_queue": _get("PMA_PBS_QUEUE", "pbs_queue"),
        "pbs_project": _get("PMA_PBS_PROJECT", "pbs_project"),
        "slurm_partition": _get("PMA_SLURM_PARTITION", "slurm_partition"),
        "slurm_account": _get("PMA_SLURM_ACCOUNT", "slurm_account"),
        "resources_scale": _get("PMA_RESOURCES_SCALE", "resources_scale"),
        "conda_env": _get("PMA_CONDA_ENV", "conda_env"),
    }

    return {
        "mode": mode,
        "hpc_env_path": "deliverables/pre_run/config/hpc.env"
        if hpc_env_path.exists() else None,
        "settings": settings,
        "s0_policy": (
            "S0 ingest runs on the login node first; on OOM/walltime (or very large "
            "inputs) it is retried as a cluster job before P2 continues."
        ),
    }


_RESOURCE_FAILURE_MARKERS = (
    "out of memory",
    "oom",
    "memoryerror",
    "cannot allocate memory",
    "killed",
    "signal 9",
    "sigkill",
    "walltime",
    "time limit",
    "exceeded memory",
    "std::bad_alloc",
    "memory allocation failed",
)


def looks_like_resource_failure(text: str) -> bool:
    """True when stderr/log output suggests OOM or walltime rather than a logic error."""
    lowered = text.lower()
    return any(marker in lowered for marker in _RESOURCE_FAILURE_MARKERS)
