"""HPC helpers — profile paths, head-job submission, scheduler detection.

Used by `executor.cli` for the `--executor` flag and the `submit` command.
Keeps cluster knowledge out of the CLI module itself.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Literal


Executor = Literal["local", "pbs", "slurm"]

# Repo root: processing-MuAgent/ — derived from this file's location.
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


def submit_head_job(
    executor: Executor,
    config_path: Path | str,
    target: str = "all",
    *,
    output_log: Path | None = None,
) -> str:
    """Submit the snakemake runner as a head-job on the chosen scheduler.

    Returns the scheduler-assigned job id (e.g. PBS "1234567.pbs" or SLURM "1234567").
    Raises CalledProcessError if submission fails.

    The head-job activates the project conda env and runs snakemake with the
    chosen profile. Per-stage child jobs are submitted by snakemake itself.
    """
    if executor not in ("pbs", "slurm"):
        raise ValueError(f"submit_head_job requires pbs|slurm; got {executor!r}")
    config_path = Path(config_path).resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"config not found: {config_path}")

    runner = RUNNER_SCRIPT[executor]
    if not runner.exists():
        raise FileNotFoundError(f"head-job script missing: {runner}")

    env_vars = {
        "PMA_CONFIG": str(config_path),
        "PMA_TARGET": target,
        "PMA_REPO_ROOT": str(REPO_ROOT),
    }

    if executor == "pbs":
        cmd = ["qsub", "-terse"]
        # Inherit the submitter's env (queue, project, notify email, etc.).
        cmd += ["-V"]
        # Plus explicit pass-through of the run-specific vars (more reliable
        # than relying on -V across all PBS Pro configurations).
        cmd += ["-v", ",".join(f"{k}={v}" for k, v in env_vars.items())]
        if output_log is not None:
            cmd += ["-o", str(output_log), "-j", "oe"]
        # Optional queue / project from env vars.
        if os.environ.get("PMA_PBS_QUEUE"):
            cmd += ["-q", os.environ["PMA_PBS_QUEUE"]]
        if os.environ.get("PMA_PBS_PROJECT"):
            cmd += ["-P", os.environ["PMA_PBS_PROJECT"]]
        cmd += [str(runner)]

    else:  # slurm
        export_list = "ALL," + ",".join(f"{k}={v}" for k, v in env_vars.items())
        cmd = ["sbatch", "--parsable", f"--export={export_list}"]
        if output_log is not None:
            cmd += ["--output", str(output_log)]
        if os.environ.get("PMA_SLURM_PARTITION"):
            cmd += ["--partition", os.environ["PMA_SLURM_PARTITION"]]
        if os.environ.get("PMA_SLURM_ACCOUNT"):
            cmd += ["--account", os.environ["PMA_SLURM_ACCOUNT"]]
        cmd += [str(runner)]

    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout.strip()


def env_diagnostics() -> dict[str, str | None]:
    """Snapshot of HPC-relevant env vars — used by `processing-muagent hpc-info`
    to show the user what's wired up.
    """
    keys = (
        "PMA_PBS_QUEUE", "PMA_PBS_PROJECT",
        "PMA_SLURM_PARTITION", "PMA_SLURM_ACCOUNT",
        "PMA_NOTIFY_EMAIL", "PMA_RESOURCES_SCALE",
        "PMA_CONDA_ENV", "PMA_LOG_DIR",
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


def write_hpc_env(path: Path | str, *, mode: Executor, settings: dict[str, str | None]) -> Path:
    """Write a source-able shell snippet with PMA_* exports for this run."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# processing-MuAgent HPC settings — source before submit/run on cluster:",
        "#   source deliverables/pre_run/config/hpc.env",
        f"# execution mode: {mode}",
        "",
    ]
    key_map = {
        "pbs_queue": "PMA_PBS_QUEUE",
        "pbs_project": "PMA_PBS_PROJECT",
        "slurm_partition": "PMA_SLURM_PARTITION",
        "slurm_account": "PMA_SLURM_ACCOUNT",
        "notify_email": "PMA_NOTIFY_EMAIL",
        "resources_scale": "PMA_RESOURCES_SCALE",
        "conda_env": "PMA_CONDA_ENV",
    }
    for field, env_key in key_map.items():
        val = settings.get(field)
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
