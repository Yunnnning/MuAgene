"""Machine-capability probing for Execution-MuAgent.

Detects what the *current* machine can do — scheduler(s), GPU presence, the conda
env manager, and the container runtime — so `configure-execution` can suggest
sensible values and `provision-env`/`validate-env` can pick the right tools. This
is the runtime-infrastructure counterpart to Processing-MuAgent's `discover_site`
(which probes scheduler queues/partitions); kept here because environment + machine
provisioning is Execution-MuAgent's responsibility.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from typing import Any


def _which(*names: str) -> str | None:
    for n in names:
        if shutil.which(n):
            return n
    return None


def _run_ok(cmd: list[str], timeout_s: int = 5) -> bool:
    try:
        return subprocess.run(cmd, capture_output=True, timeout=timeout_s).returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def detect_manager() -> str | None:
    """Preferred conda env manager on PATH (micromamba > mamba > conda)."""
    return _which("micromamba", "mamba", "conda")


def detect_container_runtime() -> str | None:
    """Container runtime on PATH (apptainer or singularity). May be gated behind a
    module load on HPC — see `environments.singularity_module` in site.config."""
    return _which("apptainer", "singularity")


def detect_scheduler() -> str:
    if shutil.which("sbatch"):
        return "slurm"
    return "local"


def gpu_present() -> bool:
    """True if this machine exposes a GPU now (nvidia-smi) or the scheduler advertises
    a gpu gres. Login nodes often have no local GPU even when the cluster does, so a
    False here is not authoritative for the cluster — only for *this* host."""
    if shutil.which("nvidia-smi") and _run_ok(["nvidia-smi", "-L"]):
        return True
    if shutil.which("sinfo"):
        try:
            out = subprocess.run(["sinfo", "-h", "-o", "%G"], capture_output=True,
                                 text=True, timeout=5).stdout
            return "gpu" in out.lower()
        except (OSError, subprocess.SubprocessError):
            return False
    return False


def probe_capabilities() -> dict[str, Any]:
    """Structured capability report for this machine."""
    return {
        "hostname": os.environ.get("HOSTNAME") or os.uname().nodename,
        "scheduler": detect_scheduler(),
        "manager": detect_manager(),
        "container_runtime": detect_container_runtime(),
        "gpu_present": gpu_present(),
    }
