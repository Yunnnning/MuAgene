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
    """Snapshot of HPC-relevant env vars — used by `processing-muagent status --hpc`
    to show the user what's wired up.
    """
    keys = (
        "PMA_PBS_QUEUE", "PMA_PBS_PROJECT",
        "PMA_SLURM_PARTITION", "PMA_SLURM_ACCOUNT",
        "PMA_NOTIFY_EMAIL", "PMA_RESOURCES_SCALE",
        "PMA_CONDA_ENV", "PMA_LOG_DIR",
    )
    return {k: os.environ.get(k) for k in keys}
