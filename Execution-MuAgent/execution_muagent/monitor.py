"""HPC job registration, monitoring, and diagnostics.

Detection and decision are separate. The watcher is cheap and deterministic: it
counts quiet check intervals and raises a flag after the stage's tolerance is
crossed. A flag is a suspicion, not a verdict. The investigator then gathers
independent evidence (CPU, memory, filesystem probe, child states, log markers)
before classifying. Kill only from confirmed dead, with a recorded reason.

Two clocks: check interval (sampling rate, same for every stage) and tolerance_n
(how many quiet intervals are allowed, derived from the stage's
progress_timeout_hint). Silence is measured in missed heartbeats, not wall-clock
minutes. A heartbeat occurs when any run-scoped file mtime advances OR the head
log grows — whichever fires first.

All scheduler calls are bounded by short subprocess timeouts so the monitor
cannot hang behind a stuck `squeue`/`qstat` call.
"""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import threading
import time
from dataclasses import asdict, dataclass, field, replace
from datetime import datetime, timezone
from enum import Enum
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

SCHEDULER_TERMINAL_STATES = SCHEDULER_FAILED_STATES | {"COMPLETED"}

class JobHealth(str, Enum):
    """Monitor state machine states for a single HPC submission."""
    HEALTHY       = "healthy"        # No stall signal; normal operation
    SUSPECT       = "suspect"        # Stall flag raised; evidence not yet gathered
    INVESTIGATING = "investigating"  # Transient: gathering evidence this iteration
    RECOVERED     = "recovered"      # Investigation found life; silence window reset
    CONFIRMED_DEAD = "confirmed_dead"  # Evidence confirmed dead; ready to kill
    FS_HANG       = "fs_hang"        # D-state / storage-degraded; killed + reported
    KILLED        = "killed"         # Cancellation sent
    DONE          = "done"           # Job exited terminal scheduler state


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
    spec_path: str | None = None
    progress_timeout_hint: float | None = None


@dataclass
class MonitorFinding:
    severity: str
    code: str
    message: str


@dataclass
class MonitorState:
    """Per-submission state carried across watch iterations."""
    health: JobHealth = JobHealth.HEALTHY
    silence_intervals: int = 0          # consecutive quiet check intervals
    tolerance_n: int = 20               # flag after this many quiet intervals (default: 20×270s=90min)
    last_progress_mtime: float | None = None  # newest progress file mtime at last check
    last_log_size: int | None = None    # head log byte-size at last check
    investigation: dict | None = None   # evidence gathered at last SUSPECT→INVESTIGATING transition
    confirmed_dead_reason: str | None = None
    previous_finding_codes: frozenset[str] | None = None
    verified_stages: frozenset[str] = field(default_factory=frozenset)  # stages whose outputs verified


@dataclass
class SiteConfig:
    """Platform description written by Processing-MuAgent, consumed by Execution-MuAgent."""
    scheduler: str
    partition: str | None = None
    account: str | None = None
    qos: str | None = None
    queue: str | None = None
    project: str | None = None
    resources_scale: int = 1
    conda_env: str | None = None
    container: str | None = None
    scratch: str | None = None


@dataclass
class StageSpec:
    """Per-stage job spec authored by Processing-MuAgent."""
    schema_version: str
    stage: str
    science_description: str
    resources: dict[str, int]
    inputs: dict[str, str]
    outputs: dict[str, str]
    progress_timeout_hint: float


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def run_monitor_dir(run_dir: Path | str) -> Path:
    return Path(run_dir) / "internal" / "hpc_monitor"


def load_site_config(path: Path | str) -> SiteConfig:
    """Load a site.config YAML file into a SiteConfig dataclass."""
    import yaml
    p = Path(path)
    data = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    scheduler = str(data.get("scheduler", "slurm"))
    sched_section: dict[str, Any] = {}
    if scheduler == "slurm":
        sched_section = data.get("slurm") or {}
    elif scheduler == "pbs":
        sched_section = data.get("pbs") or {}
    common: dict[str, Any] = data.get("common") or {}
    return SiteConfig(
        scheduler=scheduler,
        partition=sched_section.get("partition"),
        account=sched_section.get("account"),
        qos=sched_section.get("qos"),
        queue=sched_section.get("queue"),
        project=sched_section.get("project"),
        resources_scale=int(common.get("resources_scale") or 1),
        conda_env=common.get("conda_env"),
        container=common.get("container"),
        scratch=common.get("scratch"),
    )


def load_stage_spec(path: Path | str) -> StageSpec:
    """Load a per-stage job spec YAML into a StageSpec dataclass."""
    import yaml
    p = Path(path)
    data = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    return StageSpec(
        schema_version=str(data.get("schema_version", "1")),
        stage=str(data.get("stage", "")),
        science_description=str(data.get("science_description", "")),
        resources=dict(data.get("resources") or {}),
        inputs=dict(data.get("inputs") or {}),
        outputs=dict(data.get("outputs") or {}),
        progress_timeout_hint=float(data.get("progress_timeout_hint", 90)),
    )


_HDF5_SIGNATURE = b"\x89HDF\r\n\x1a\n"
# HDF5 superblock may start at offset 0 or a later power-of-two offset.
_HDF5_SIGNATURE_OFFSETS = (0, 512, 1024, 2048, 4096, 8192)


def _looks_like_hdf5(path: Path) -> bool:
    """Dependency-free HDF5 check: the 8-byte signature at a superblock offset."""
    try:
        with path.open("rb") as fh:
            for off in _HDF5_SIGNATURE_OFFSETS:
                fh.seek(off)
                if fh.read(8) == _HDF5_SIGNATURE:
                    return True
    except OSError:
        return False
    return False


def _looks_like_parquet(path: Path) -> bool:
    """Dependency-free parquet check: the PAR1 magic at head and tail."""
    try:
        size = path.stat().st_size
        if size < 8:
            return False
        with path.open("rb") as fh:
            head = fh.read(4)
            fh.seek(-4, os.SEEK_END)
            tail = fh.read(4)
        return head == b"PAR1" and tail == b"PAR1"
    except OSError:
        return False


def verify_output_file(path: Path | str) -> tuple[bool, str]:
    """Verify that an output artifact is present and a complete, loadable file.

    Lightweight integrity check (not a full in-memory load) suitable for running
    on every monitor interval. Real loaders (h5py / pyarrow) are used when
    available in the runtime env; otherwise we fall back to dependency-free
    structural checks. Because stages write outputs atomically (write to /tmp,
    fsync, then os.rename), a file present at its final path is complete — so a
    valid signature plus non-zero size is a strong correctness signal.

    Returns (ok, reason). reason is a short machine token, e.g. "valid_hdf5",
    "missing", "empty", "corrupt_hdf5".
    """
    p = Path(path)
    if not p.exists():
        return False, "missing"
    try:
        if p.stat().st_size == 0:
            return False, "empty"
    except OSError as exc:
        return False, f"stat_error:{exc}"

    suffix = p.suffix.lower()
    if suffix in (".h5ad", ".h5mu", ".h5"):
        try:
            import h5py  # type: ignore[import-not-found]

            with h5py.File(str(p), "r") as fh:
                keys = set(fh.keys())
            if suffix == ".h5mu":
                ok = "mod" in keys
            elif suffix == ".h5ad":
                ok = bool({"X", "obs", "var"} & keys)
            else:
                ok = bool(keys)
            return (True, "valid_hdf5") if ok else (False, "hdf5_missing_groups")
        except ImportError:
            return (True, "hdf5_signature") if _looks_like_hdf5(p) else (False, "corrupt_hdf5")
        except Exception as exc:  # truncated / unreadable HDF5
            return False, f"corrupt_hdf5:{type(exc).__name__}"
    if suffix == ".parquet":
        try:
            import pyarrow.parquet as pq  # type: ignore[import-not-found]

            pq.read_metadata(str(p))
            return True, "valid_parquet"
        except ImportError:
            return (True, "parquet_magic") if _looks_like_parquet(p) else (False, "corrupt_parquet")
        except Exception as exc:
            return False, f"corrupt_parquet:{type(exc).__name__}"
    if suffix == ".json":
        try:
            json.loads(p.read_text(encoding="utf-8"))
            return True, "valid_json"
        except Exception as exc:
            return False, f"corrupt_json:{type(exc).__name__}"
    # Text sentinels and any other type: non-zero size already confirmed above.
    return True, "nonempty"


def verify_stage_outputs(submission: Submission) -> dict[str, dict[str, tuple[bool, str]]]:
    """Verify declared outputs for every per-stage spec under internal/stage_meta/.

    Returns {stage: {output_name: (ok, reason)}}. Only declared outputs are
    checked; the head_job spec (no outputs) is skipped.
    """
    stage_meta_dir = Path(submission.run_dir) / "internal" / "stage_meta"
    results: dict[str, dict[str, tuple[bool, str]]] = {}
    if not stage_meta_dir.exists():
        return results
    for meta_path in sorted(stage_meta_dir.glob("*.yaml")):
        if meta_path.stem == "head_job":
            continue
        try:
            spec = load_stage_spec(meta_path)
        except Exception:
            continue
        if not spec.outputs:
            continue
        results[spec.stage] = {
            name: verify_output_file(path_str) for name, path_str in spec.outputs.items()
        }
    return results


def validate_spec(spec: StageSpec, site_config: SiteConfig) -> list[str]:
    """Return a list of validation error strings (empty = spec is workable)."""
    errors: list[str] = []
    if not spec.stage:
        errors.append("spec.stage is empty")
    if spec.resources.get("cpus", 0) <= 0:
        errors.append(f"spec.resources.cpus must be > 0, got {spec.resources.get('cpus')}")
    if spec.resources.get("mem_mb", 0) <= 0:
        errors.append(f"spec.resources.mem_mb must be > 0, got {spec.resources.get('mem_mb')}")
    if spec.resources.get("walltime_min", 0) <= 0:
        errors.append(f"spec.resources.walltime_min must be > 0, got {spec.resources.get('walltime_min')}")
    if site_config.scheduler not in ("slurm", "pbs"):
        errors.append(f"site_config.scheduler must be slurm or pbs, got {site_config.scheduler!r}")
    for key, path in spec.inputs.items():
        if path and not Path(path).exists():
            errors.append(f"spec input {key!r} not found: {path}")
    return errors


def render_submission_script(
    spec: StageSpec,
    site_config: SiteConfig,
    repo_root: Path | str,
    run_dir: Path | str,
    log_path: Path | str,
    target: str,
) -> str:
    """Render a scheduler submission script from a stage spec + site.config.

    Maps resource hints to scheduler directives, adds partition/account/QOS,
    wraps the command in a container invocation when site_config.container is set,
    and invokes launch_runner.sh with the correct environment.
    """
    repo_root = Path(repo_root)
    run_dir_s = str(Path(run_dir).resolve())
    log_path_s = str(Path(log_path).resolve())
    walltime_min = spec.resources.get("walltime_min", 60)
    cpus = spec.resources.get("cpus", 1)
    mem_mb = spec.resources.get("mem_mb", 4000)
    conda_env = site_config.conda_env or "grn"

    if site_config.scheduler == "slurm":
        hh, mm = divmod(walltime_min, 60)
        lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name=pma_{spec.stage}",
            f"#SBATCH --output={log_path_s}",
            f"#SBATCH --cpus-per-task={cpus}",
            f"#SBATCH --mem={mem_mb}M",
            f"#SBATCH --time={hh:02d}:{mm:02d}:00",
        ]
        if site_config.partition:
            lines.append(f"#SBATCH --partition={site_config.partition}")
        if site_config.account:
            lines.append(f"#SBATCH --account={site_config.account}")
        if site_config.qos:
            lines.append(f"#SBATCH --qos={site_config.qos}")
    else:  # pbs
        hh, mm = divmod(walltime_min, 60)
        lines = [
            "#!/bin/bash",
            f"#PBS -N pma_{spec.stage}",
            f"#PBS -o {log_path_s}",
            "#PBS -j oe",
            f"#PBS -l select=1:ncpus={cpus}:mem={mem_mb}mb",
            f"#PBS -l walltime={hh:02d}:{mm:02d}:00",
        ]
        if site_config.queue:
            lines.append(f"#PBS -q {site_config.queue}")
        if site_config.project:
            lines.append(f"#PBS -P {site_config.project}")

    run_dir_path = Path(run_dir).resolve()
    lines += [
        "",
        f"export PMA_CONFIG={run_dir_path / 'deliverables' / 'pre_run' / 'config' / 'run.yaml'}",
        f"export PMA_TARGET={target}",
        f"export PMA_REPO_ROOT={repo_root}",
    ]
    # Export scheduler-specific vars so launch_runner.sh adds them as
    # --default-resources when it detects the cluster profile.
    if conda_env:
        lines.append(f"export PMA_CONDA_ENV={conda_env}")
    if site_config.scheduler == "slurm":
        if site_config.partition:
            lines.append(f"export PMA_SLURM_PARTITION={site_config.partition}")
        if site_config.account:
            lines.append(f"export PMA_SLURM_ACCOUNT={site_config.account}")
    elif site_config.scheduler == "pbs":
        if site_config.queue:
            lines.append(f"export PMA_PBS_QUEUE={site_config.queue}")
        if site_config.project:
            lines.append(f"export PMA_PBS_PROJECT={site_config.project}")
    lines.append("")

    profile = repo_root / "workflow" / "profiles" / site_config.scheduler
    launch = str(repo_root / "scripts" / "launch_runner.sh")
    launch_args = f"--configfile $PMA_CONFIG --profile {profile} --jobs 8 $PMA_TARGET"
    if site_config.container:
        lines.append(
            f"apptainer exec --bind {run_dir_s}:{run_dir_s} "
            f"{site_config.container} bash {launch} {launch_args}"
        )
    else:
        lines += [
            f"source $(conda info --base)/etc/profile.d/conda.sh",
            f"conda activate {conda_env}",
            f"bash {launch} {launch_args}",
        ]

    return "\n".join(lines) + "\n"


def submit_from_spec(
    spec: StageSpec,
    site_config: SiteConfig,
    run_dir: Path | str,
    repo_root: Path | str,
    log_path: Path | str,
    target: str,
    timeout_s: int = 60,
) -> dict[str, Any]:
    """Render a submission script, write it to disk, and submit via sbatch/qsub.

    Returns a dict with job_id, script_path, rejected_as ('policy'|'transient'|None),
    stdout, and stderr. The caller registers the submission and starts monitoring.
    """
    run_dir = Path(run_dir).resolve()
    scripts_dir = run_monitor_dir(run_dir) / "scripts"
    scripts_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    script_path = scripts_dir / f"{spec.stage}_{stamp}.sh"

    script_text = render_submission_script(spec, site_config, repo_root, run_dir, log_path, target)
    script_path.write_text(script_text, encoding="utf-8")
    script_path.chmod(0o755)

    _POLICY_MARKERS = (
        "invalid partition", "invalid account", "invalid qos",
        "time limit", "exceeds", "not found", "error in batch script",
        "invalid node", "unrecognized", "no partition specified",
    )

    def _try_submit() -> tuple[int | None, str, str]:
        if site_config.scheduler == "slurm":
            cmd = ["sbatch", "--parsable", str(script_path)]
        else:
            cmd = ["qsub", "-terse", str(script_path)]
        return _run_cmd(cmd, timeout_s)

    rc, out, err = _try_submit()
    rejected_as: str | None = None

    if rc != 0:
        combined = (out + err).lower()
        if any(m in combined for m in _POLICY_MARKERS):
            rejected_as = "policy"
        else:
            # Transient failure — retry up to 2 times with a short backoff.
            for _ in range(2):
                time.sleep(10)
                rc, out, err = _try_submit()
                if rc == 0:
                    break
            if rc != 0:
                rejected_as = "transient"

    job_id = out.strip().split()[0] if rc == 0 and out.strip() else ""
    return {
        "job_id": job_id,
        "script_path": str(script_path),
        "returncode": rc,
        "stdout": out,
        "stderr": err,
        "rejected_as": rejected_as,
    }


def append_execution_manifest(run_dir: Path | str, entry: dict[str, Any]) -> None:
    """Append one entry to the execution manifest (job_id, spec, script, outputs)."""
    manifest = run_monitor_dir(run_dir) / "execution_manifest.jsonl"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry, sort_keys=True) + "\n")


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
    """Return run-local log roots first, with legacy repo-level roots as fallback.

    Repo-level ``logs/`` is intentionally excluded: concurrent runs share that
    directory and it poisons stale-progress detection.
    """
    run_dir = Path(submission.run_dir)
    repo_root = Path(submission.repo_root)
    return [
        run_dir / "internal" / "snakemake" / ".snakemake" / "log",
        run_dir / "internal" / "snakemake" / ".snakemake" / "slurm_logs",
        repo_root / ".snakemake" / "log",
        repo_root / ".snakemake" / "slurm_logs",
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


def _is_run_scoped_progress_path(submission: Submission, path: Path) -> bool:
    run_dir = Path(submission.run_dir).resolve()
    try:
        path.resolve().relative_to(run_dir)
        return True
    except ValueError:
        return path.resolve() == Path(submission.log_path).resolve()


def discover_progress_files(submission: Submission) -> list[Path]:
    """Progress files used for stale-run detection (scoped to this run)."""
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
    scoped = [p for p in files if _is_run_scoped_progress_path(submission, p)]
    return sorted(set(scoped), key=lambda p: str(p))


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
    submitted_after = max(0.0, _submitted_epoch(submission) - 60.0)
    ids: set[str] = set()
    slurm_logs = (
        Path(submission.run_dir) / "internal" / "snakemake" / ".snakemake" / "slurm_logs"
    )
    if slurm_logs.exists():
        for path in slurm_logs.rglob("*.log"):
            if not path.is_file() or path.stat().st_mtime < submitted_after:
                continue
            if path.stem.isdigit():
                norm = _normalize_job_id(path.stem)
                if norm and norm != head_id:
                    ids.add(norm)
    for path in discover_log_files(submission):
        for job_id in parse_job_ids_from_log(path):
            norm = _normalize_job_id(job_id)
            if norm and norm != head_id:
                ids.add(norm)
    return sorted(ids)


def _child_slurm_log(submission: Submission, job_id: str) -> Path | None:
    root = Path(submission.run_dir) / "internal" / "snakemake" / ".snakemake" / "slurm_logs"
    if not root.exists():
        return None
    matches = [p for p in root.rglob(f"{_normalize_job_id(job_id)}.log") if p.is_file()]
    return _newest_file(matches)


def _log_indicates_storage_hang(path: Path) -> bool:
    """True if a Snakemake child log shows finished but stuck at output storage."""
    if not path.exists() or not path.is_file():
        return False
    tail = _read_tail(path, 4000)
    return "Finished jobid:" in tail and "Storing output in storage." in tail


def _child_storage_hang_ids(
    submission: Submission,
    child_job_ids: list[str],
    scheduler_timeout_s: int,
) -> list[str]:
    """Return child job IDs that are RUNNING but stuck at Snakemake output storage."""
    hung: list[str] = []
    for child_id in child_job_ids:
        sched = query_scheduler(submission.executor, child_id, scheduler_timeout_s)
        state = str(sched.get("state") or "").split()[0]
        if state not in SCHEDULER_RUNNING_STATES:
            continue
        child_log = _child_slurm_log(submission, child_id)
        if child_log is not None and _log_indicates_storage_hang(child_log):
            hung.append(child_id)
    return hung


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


def evaluate_snapshot_definitive(
    snapshot: dict[str, Any],
) -> list[MonitorFinding]:
    """Return only findings whose meaning is unambiguous regardless of elapsed time.

    Stall detection (silence counting) is handled by the state machine in
    monitor_once; this function only raises findings that warrant immediate
    attention independent of any silence window.
    """
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
    if not (snapshot.get("latest_progress_file") or {}).get("exists") and state in SCHEDULER_RUNNING_STATES:
        findings.append(MonitorFinding("warning", "no_progress_files", "No run-scoped progress files found yet."))
    if state == "COMPLETING":
        findings.append(
            MonitorFinding(
                "warning",
                "scheduler_completing",
                "Job is in COMPLETING; on this cluster that can indicate epilog or NFS cleanup stall.",
            )
        )
    return findings


def _finding_codes(findings: list[MonitorFinding]) -> frozenset[str]:
    return frozenset(f.code for f in findings)


def _scheduler_state(info: dict[str, str | None]) -> str:
    return str(info.get("state") or "").split()[0]


def submission_jobs_active(
    submission: Submission,
    snapshot: dict[str, Any],
    scheduler_timeout_s: int = 5,
) -> bool:
    """True if the registered head job or any discovered child is still active."""
    head_state = _scheduler_state(
        query_scheduler(submission.executor, submission.job_id, scheduler_timeout_s)
    )
    if head_state in SCHEDULER_RUNNING_STATES or head_state == "PENDING":
        return True
    for child_id in snapshot.get("child_job_ids") or []:
        child_state = _scheduler_state(
            query_scheduler(submission.executor, child_id, scheduler_timeout_s)
        )
        if child_state in SCHEDULER_RUNNING_STATES or child_state == "PENDING":
            return True
    return False


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


def render_report(
    snapshot: dict[str, Any],
    findings: list[MonitorFinding],
    cancel_result: dict[str, Any] | None,
    monitor_state: "MonitorState | None" = None,
) -> str:
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
        "## Monitor State",
        "",
    ]
    if monitor_state is not None:
        lines.append(f"- Health: `{monitor_state.health.value}`")
        lines.append(f"- Silence: `{monitor_state.silence_intervals} / {monitor_state.tolerance_n} intervals`")
        verified = sorted(monitor_state.verified_stages)
        lines.append(f"- Verified stages: `{', '.join(verified) if verified else 'none'}`")
        if monitor_state.confirmed_dead_reason:
            lines.append(f"- Confirmed-dead reason: `{monitor_state.confirmed_dead_reason}`")
        if monitor_state.investigation:
            inv = monitor_state.investigation
            cpu = inv.get("compute_cpu_s")
            rss = inv.get("memory_rss_mb")
            fs = inv.get("filesystem_responsive")
            lines.append(
                f"- Investigation: filesystem_responsive={fs}, "
                f"compute_cpu_s={cpu}, memory_rss_mb={rss}, "
                f"child_storage_hang_ids={inv.get('child_storage_hang_ids')}"
            )
    else:
        lines.append("- No monitor state recorded.")
    lines.extend(["", "## Findings", ""])
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


def _monitor_state_dict(monitor_state: "MonitorState") -> dict[str, Any]:
    return {
        "health": monitor_state.health.value,
        "silence_intervals": monitor_state.silence_intervals,
        "tolerance_n": monitor_state.tolerance_n,
        "investigation": monitor_state.investigation,
        "confirmed_dead_reason": monitor_state.confirmed_dead_reason,
        "verified_stages": sorted(monitor_state.verified_stages),
    }


def _persist_snapshot(
    submission: Submission,
    snapshot: dict[str, Any],
    monitor_state: "MonitorState | None",
) -> Path:
    """Write latest_snapshot.json (snapshot + monitor_state) for Processing to read.

    Called on every check — including healthy, finding-free ones — so Processing's
    hpc-status never reads stale health.
    """
    full_state: dict[str, Any] = dict(snapshot)
    if monitor_state is not None:
        full_state["monitor_state"] = _monitor_state_dict(monitor_state)
    state_path = run_monitor_dir(submission.run_dir) / "latest_snapshot.json"
    state_path.parent.mkdir(parents=True, exist_ok=True)
    state_path.write_text(json.dumps(full_state, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return state_path


def write_report(
    submission: Submission,
    snapshot: dict[str, Any],
    findings: list[MonitorFinding],
    cancel_result: dict[str, Any] | None = None,
    *,
    monitor_state: "MonitorState | None" = None,
    record_history: bool = True,
) -> Path:
    report_dir = run_monitor_dir(submission.run_dir) / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    text = render_report(snapshot, findings, cancel_result, monitor_state)
    latest = run_monitor_dir(submission.run_dir) / "latest_report.md"
    latest.write_text(text, encoding="utf-8")
    # Persist full snapshot + monitor state so Processing-MuAgent can read health.
    _persist_snapshot(submission, snapshot, monitor_state)
    if not record_history:
        return latest
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    report_path = report_dir / f"{submission.job_id}_{stamp}.md"
    report_path.write_text(text, encoding="utf-8")
    return report_path


def _resolve_stale_minutes(submission: Submission, stale_minutes: float) -> float:
    """Return stale_minutes, preferring the spec's progress_timeout_hint when set."""
    if submission.progress_timeout_hint is not None:
        return float(submission.progress_timeout_hint)
    if submission.spec_path:
        try:
            spec = load_stage_spec(submission.spec_path)
            return spec.progress_timeout_hint
        except Exception:
            pass
    return stale_minutes


def _resolve_tolerance_n(submission: Submission, interval_s: float, fallback_minutes: float = 90.0) -> int:
    """Return how many consecutive quiet intervals constitute a stall.

    Converts the stage's progress_timeout_hint (minutes) to an interval count.
    The fallback is 90 min — conservative enough that an absent hint never causes
    a premature kill.
    """
    stale_minutes = _resolve_stale_minutes(submission, fallback_minutes)
    return max(1, math.ceil(stale_minutes * 60.0 / max(1.0, interval_s)))


def _watcher_update(
    snapshot: dict[str, Any],
    state: MonitorState,
) -> tuple[MonitorState, bool]:
    """Update silence counter from the latest snapshot.

    A heartbeat fires when any run-scoped file mtime advanced OR the head log
    grew since the previous check. On a heartbeat silence resets to 0; otherwise
    it increments. Returns (new_state, flag_raised) where flag_raised is True
    once silence_intervals reaches tolerance_n.
    """
    latest = snapshot.get("latest_progress_file") or {}
    raw_mtime = latest.get("mtime")
    current_mtime: float | None = float(raw_mtime) if raw_mtime is not None else None

    head_log = snapshot.get("head_log") or {}
    raw_size = head_log.get("size")
    current_log_size: int | None = int(raw_size) if raw_size is not None else None

    file_heartbeat = (
        current_mtime is not None
        and (state.last_progress_mtime is None or current_mtime > state.last_progress_mtime)
    )
    log_heartbeat = (
        current_log_size is not None
        and (state.last_log_size is None or current_log_size > state.last_log_size)
    )
    heartbeat = file_heartbeat or log_heartbeat

    new_silence = 0 if heartbeat else state.silence_intervals + 1
    flag_raised = new_silence >= state.tolerance_n

    new_state = replace(
        state,
        silence_intervals=new_silence,
        last_progress_mtime=current_mtime if current_mtime is not None else state.last_progress_mtime,
        last_log_size=current_log_size if current_log_size is not None else state.last_log_size,
    )
    return new_state, flag_raised


def _probe_sstat(executor: str, job_id: str, timeout_s: int = 5) -> dict[str, Any]:
    """Query live job CPU/memory via sstat (SLURM only). Best-effort."""
    if executor != "slurm":
        return {"error": "not_slurm"}
    if not shutil.which("sstat"):
        return {"error": "sstat_not_found"}
    jid = _normalize_job_id(job_id)
    rc, out, err = _run_cmd(
        ["sstat", "-j", f"{jid}.batch", "--format=AveCPU,MaxRSS,JobID", "-n", "-P"],
        timeout_s,
    )
    if rc != 0 or not out.strip():
        return {"error": err or f"sstat rc={rc}"}
    parts = out.splitlines()[0].split("|")
    ave_cpu_raw = parts[0].strip() if len(parts) > 0 else ""
    max_rss_raw = parts[1].strip() if len(parts) > 1 else ""

    ave_cpu_s: float | None = None
    if ave_cpu_raw:
        try:
            fields = ave_cpu_raw.split(":")
            if len(fields) == 3:
                ave_cpu_s = int(fields[0]) * 3600 + int(fields[1]) * 60 + float(fields[2])
        except (ValueError, IndexError):
            pass

    max_rss_mb: float | None = None
    if max_rss_raw:
        try:
            val = max_rss_raw.upper()
            if val.endswith("K"):
                max_rss_mb = float(val[:-1]) / 1024.0
            elif val.endswith("M"):
                max_rss_mb = float(val[:-1])
            elif val.endswith("G"):
                max_rss_mb = float(val[:-1]) * 1024.0
            else:
                max_rss_mb = float(val) / (1024.0 * 1024.0)  # assume bytes
        except ValueError:
            pass

    return {"ave_cpu_s": ave_cpu_s, "max_rss_mb": max_rss_mb}


def _probe_filesystem(path: str | Path, timeout_s: float = 5.0) -> bool:
    """True if os.stat(path) returns within timeout_s (filesystem is responsive)."""
    result = [False]

    def _check() -> None:
        try:
            os.stat(str(path))
            result[0] = True
        except OSError:
            result[0] = True  # stat responded even with error

    t = threading.Thread(target=_check, daemon=True)
    t.start()
    t.join(timeout_s)
    return result[0]


def investigate_suspect(
    submission: Submission,
    snapshot: dict[str, Any],
    *,
    scheduler_timeout_s: int = 5,
) -> dict[str, Any]:
    """Gather independent evidence about a stall suspect.

    Called once when the watcher transitions the job to SUSPECT. Evidence is
    collected synchronously (scheduler query already in snapshot; sstat and
    filesystem probe run fresh). The caller (classify_investigation) applies
    rules to this evidence to conclude recovered, confirmed_dead, or fs_hang.
    """
    scheduler = snapshot.get("scheduler") or {}
    scheduler_state = str(scheduler.get("state") or "").split()[0]
    error_markers = list(snapshot.get("error_markers") or [])

    sstat = _probe_sstat(submission.executor, submission.job_id, scheduler_timeout_s)
    fs_responsive = _probe_filesystem(submission.run_dir, timeout_s=5.0)

    child_ids = list(snapshot.get("child_job_ids") or [])
    hung_children = _child_storage_hang_ids(submission, child_ids, scheduler_timeout_s)

    return {
        "scheduler_state": scheduler_state,
        "error_markers": error_markers,
        "compute_cpu_s": sstat.get("ave_cpu_s"),
        "memory_rss_mb": sstat.get("max_rss_mb"),
        "sstat_error": sstat.get("error"),
        "filesystem_responsive": fs_responsive,
        "child_storage_hang_ids": hung_children,
    }


def classify_investigation(evidence: dict[str, Any]) -> tuple[str, str]:
    """Classify stall evidence. Returns (verdict, reason).

    verdict is one of: "confirmed_dead", "recovered", "fs_hang".

    Rules applied in order; first match wins. Asymmetry is intentional:
    positive evidence of life (CPU active, memory large) overrides silence;
    confirmed_dead requires all-silent signals plus a responsive filesystem
    (ruling out D-state). When evidence is incomplete, recover rather than kill.
    """
    scheduler_state = str(evidence.get("scheduler_state") or "")
    error_markers = list(evidence.get("error_markers") or [])
    hung_children = list(evidence.get("child_storage_hang_ids") or [])
    fs_responsive = evidence.get("filesystem_responsive")
    cpu_s = evidence.get("compute_cpu_s")
    rss_mb = evidence.get("memory_rss_mb")

    # 1. Scheduler reports definitive failure
    if scheduler_state in SCHEDULER_FAILED_STATES:
        return "confirmed_dead", f"scheduler_state={scheduler_state}"

    # 2. Error markers in logs — unambiguous workflow failure
    if error_markers:
        return "confirmed_dead", f"workflow_error_markers={','.join(error_markers)}"

    # 3. Child stuck at output storage — filesystem hang pattern
    if hung_children:
        return "fs_hang", f"child_storage_hang child_ids={','.join(hung_children)}"

    # 4. Filesystem probe timed out — D-state / degraded storage
    if fs_responsive is False:
        return "fs_hang", "filesystem_probe_timed_out"

    # 5. CPU activity reported by sstat — job is computing
    if cpu_s is not None and cpu_s > 1.0:
        return "recovered", f"compute_active cpu_s={cpu_s:.1f}"

    # 6. Memory still large — live process proxy
    if rss_mb is not None and rss_mb > 100.0:
        return "recovered", f"memory_active rss_mb={rss_mb:.0f}"

    # 7. Filesystem responsive + scheduler RUNNING + no CPU/memory evidence
    # All silence signals confirmed with a working filesystem → dead
    if fs_responsive is True and scheduler_state in SCHEDULER_RUNNING_STATES:
        return "confirmed_dead", "no_activity_with_responsive_filesystem"

    # 8. Inconclusive evidence — conservative: assume alive
    return "recovered", "insufficient_evidence_assume_alive"


def monitor_once(
    submission: Submission,
    state: MonitorState,
    *,
    scheduler_timeout_s: int = 5,
    kill_on_hang: bool = True,
) -> tuple[dict[str, Any], list[MonitorFinding], dict[str, Any] | None, Path | None, MonitorState]:
    """Single monitoring check driving the state machine.

    States HEALTHY / RECOVERED run the watcher (cheap: file mtime + log size).
    SUSPECT triggers investigation and classification. An unhealthy verdict
    (CONFIRMED_DEAD or FS_HANG) kills and reports — Execution never holds for a
    human and never resubmits; Processing-MuAgent owns recovery.

    Each check also verifies any per-stage outputs that have appeared and emits a
    one-time ``stage_output_verified`` progress finding, so Processing sees both
    normal progress and unhealthy flags.
    """
    snapshot = collect_snapshot(submission, scheduler_timeout_s)
    findings: list[MonitorFinding] = []
    cancel_result: dict[str, Any] | None = None

    # --- Per-step output verification (normal-progress reporting) ---
    newly_verified: list[str] = []
    for stage, outputs in verify_stage_outputs(submission).items():
        if stage in state.verified_stages:
            continue
        if outputs and all(ok for ok, _ in outputs.values()):
            newly_verified.append(stage)
    if newly_verified:
        state = replace(state, verified_stages=state.verified_stages | frozenset(newly_verified))
        for stage in sorted(newly_verified):
            findings.append(MonitorFinding(
                "info", "stage_output_verified",
                f"Stage {stage!r} outputs verified complete and loadable.",
            ))

    # --- Definitive signals (always checked, bypass silence state machine) ---
    scheduler = snapshot.get("scheduler") or {}
    sched_state = str(scheduler.get("state") or "").split()[0]

    definitive = evaluate_snapshot_definitive(snapshot)
    # Terminal definitive findings immediately override health
    for f in definitive:
        if f.code in ("scheduler_failed", "workflow_error_marker"):
            if state.health not in (JobHealth.KILLED, JobHealth.DONE):
                state = replace(state, health=JobHealth.CONFIRMED_DEAD,
                                confirmed_dead_reason=f"{f.code}: {f.message}")
    findings.extend(definitive)

    # --- State machine ---
    if state.health in (JobHealth.HEALTHY, JobHealth.RECOVERED):
        # Only run watcher when there is at least one progress file to track
        has_progress = bool((snapshot.get("latest_progress_file") or {}).get("exists"))
        if has_progress:
            state, flag_raised = _watcher_update(snapshot, state)
            if flag_raised:
                state = replace(state, health=JobHealth.SUSPECT)
                findings.append(MonitorFinding(
                    "warning", "stall_suspected",
                    f"No progress for {state.silence_intervals} consecutive check intervals "
                    f"(tolerance: {state.tolerance_n}). Entering investigation.",
                ))

    elif state.health == JobHealth.SUSPECT:
        state = replace(state, health=JobHealth.INVESTIGATING)
        evidence = investigate_suspect(submission, snapshot, scheduler_timeout_s=scheduler_timeout_s)
        state = replace(state, investigation=evidence)
        verdict, reason = classify_investigation(evidence)

        if verdict == "confirmed_dead":
            state = replace(state, health=JobHealth.CONFIRMED_DEAD, confirmed_dead_reason=reason)
            findings.append(MonitorFinding("error", "stall_confirmed",
                f"Investigation concluded confirmed dead: {reason}"))
        elif verdict == "fs_hang":
            state = replace(state, health=JobHealth.FS_HANG)
            findings.append(MonitorFinding("error", "filesystem_hang_suspected",
                f"Filesystem-related hang detected: {reason}. "
                "Killing and reporting to Processing-MuAgent for recovery."))
        else:  # recovered
            state = replace(state, health=JobHealth.HEALTHY, silence_intervals=0)
            findings.append(MonitorFinding("warning", "stall_recovered",
                f"Investigation found evidence of life: {reason}. Resuming monitoring."))

    # INVESTIGATING is transient (only set and resolved within a single call above)

    # --- Kill path ---
    # An unhealthy verdict (confirmed dead OR filesystem hang) is killed for
    # cleanup and reported. Execution never holds for a human and never
    # resubmits — Processing-MuAgent reads the report, escalates to the human,
    # fixes, and resubmits.
    child_job_ids = list(snapshot.get("child_job_ids") or [])

    if state.health in (JobHealth.CONFIRMED_DEAD, JobHealth.FS_HANG) and kill_on_hang:
        cancel_result = cancel_submission_jobs(submission, child_job_ids, scheduler_timeout_s)
        cancel_result["confirmed_dead_reason"] = (
            state.confirmed_dead_reason if state.health == JobHealth.CONFIRMED_DEAD else "filesystem_hang"
        )
        state = replace(state, health=JobHealth.KILLED)

    # --- Report ---
    # latest_snapshot.json is always refreshed (normal + unhealthy) so Processing's
    # hpc-status never reads stale state. latest_report.md / reports/ history is
    # written only when there are findings or a kill action.
    finding_codes = _finding_codes(findings)
    report_path: Path | None = None
    if findings or cancel_result:
        record_history = (
            cancel_result is not None
            or state.previous_finding_codes is None
            or finding_codes != state.previous_finding_codes
        )
        report_path = write_report(
            submission,
            snapshot,
            findings,
            cancel_result,
            monitor_state=state,
            record_history=record_history,
        )
    else:
        _persist_snapshot(submission, snapshot, state)

    new_state = replace(state, previous_finding_codes=finding_codes)
    return snapshot, findings, cancel_result, report_path, new_state


def validate_terminal_outputs(submission: Submission) -> list[MonitorFinding]:
    """After a head-job reaches a terminal COMPLETED state, verify per-stage outputs.

    Reads every stage metadata YAML in internal/stage_meta/ and properly verifies
    each declared output via verify_output_file (a loadable, non-truncated file —
    not merely non-empty). Returns an error finding for each missing, empty, or
    corrupt artifact. An empty return means all outputs verified.
    """
    findings: list[MonitorFinding] = []
    for stage, outputs in verify_stage_outputs(submission).items():
        for name, (ok, reason) in outputs.items():
            if not ok:
                findings.append(MonitorFinding(
                    severity="error",
                    code="output_missing",
                    message=(
                        f"Stage {stage!r} expected output {name!r} failed "
                        f"verification ({reason})."
                    ),
                ))
    return findings


def monitor_watch(
    submission: Submission,
    *,
    interval_s: float = 270.0,
    stale_minutes: float = 90.0,
    scheduler_timeout_s: int = 5,
    kill_on_hang: bool = True,
    max_checks: int | None = None,
) -> Path | None:
    """Poll until all jobs exit, driving the state machine across iterations."""
    tolerance_n = _resolve_tolerance_n(submission, interval_s, stale_minutes)
    state = MonitorState(tolerance_n=tolerance_n)

    checks = 0
    last_report: Path | None = None

    while True:
        snapshot, findings, cancel_result, report_path, state = monitor_once(
            submission,
            state,
            scheduler_timeout_s=scheduler_timeout_s,
            kill_on_hang=kill_on_hang,
        )
        if report_path is not None:
            last_report = report_path

        if not submission_jobs_active(submission, snapshot, scheduler_timeout_s):
            break

        checks += 1
        if max_checks is not None and checks >= max_checks:
            break
        time.sleep(max(5.0, interval_s))

    # D6: on clean terminal exit validate all per-stage output artifacts.
    if last_report is None:
        val_findings = validate_terminal_outputs(submission)
        if val_findings:
            last_report = write_report(
                submission,
                {},
                val_findings,
                None,
                monitor_state=state,
                record_history=True,
            )
    return last_report


def parse_job_ids_from_log(path: Path | str) -> list[str]:
    """Extract scheduler ids from a Snakemake/head-job log."""
    p = Path(path)
    ids: set[str] = set()
    if p.suffix == ".log" and p.stem.isdigit():
        ids.add(p.stem)
    if not p.exists():
        return sorted(ids)
    text = p.read_text(encoding="utf-8", errors="replace")
    patterns = [
        r"external jobid '([0-9]+)'",
        r"SLURM jobid ([0-9]+)",
        r"Submitted batch job ([0-9]+)",
        r"Submitted (?:slurm )?head-job: ([0-9]+)",
    ]
    for pat in patterns:
        ids.update(re.findall(pat, text))
    return sorted(ids)
