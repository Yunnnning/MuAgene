"""Compute-device resolution + GPU capability gate for stage code.

The operator selects ``compute.device`` (cpu|gpu) once; the submit path exports
``PMA_DEVICE`` onto each GPU-capable child job. Which stages are GPU-capable is
declared in ONE place — ``_GPU_CAPABLE`` in workflow/resources.smk (which also
documents the two-edit contract for adding a stage); this module is stage-agnostic
and keeps no second list. Stage code calls :func:`use_gpu` to branch a heavy op to
its rapids-singlecell drop-in, then records the device actually used in provenance.

Never silently degrade (see the team's execution-vs-scientific-errors rule): if gpu
is requested but unavailable at runtime, raise loudly so it is fixed at the env/
allocation level — unless ``PMA_DEVICE_FALLBACK=1`` explicitly permits a CPU fallback
(then warn + record). The preflight (Execution-MuAgent validate-env) normally catches
a missing GPU env before submit; this is the in-job last line of defence.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
from typing import Any


class GpuUnavailableError(RuntimeError):
    """Raised when compute.device=gpu is requested but no usable GPU/env is present."""


def requested_device() -> str:
    """The operator's requested device for this job ('cpu' | 'gpu')."""
    return (os.environ.get("PMA_DEVICE", "cpu") or "cpu").strip().lower()


def _rapids_available() -> bool:
    return importlib.util.find_spec("rapids_singlecell") is not None


def _gpu_visible() -> bool:
    try:
        import cupy  # type: ignore[import-not-found]

        return cupy.cuda.runtime.getDeviceCount() > 0
    except Exception:
        return False


def gpu_capable() -> tuple[bool, str]:
    """(is_usable, reason). Usable means rapids-singlecell imports AND a CUDA device
    is visible — both are required to actually run a GPU op."""
    if not _rapids_available():
        return False, "rapids_singlecell not importable in this environment"
    if not _gpu_visible():
        return False, "no CUDA device visible (cupy.cuda.runtime.getDeviceCount()==0)"
    return True, "rapids-singlecell + CUDA device available"


def use_gpu(*, run_dir: str | Path | None = None, stage: str | None = None) -> bool:
    """Whether this stage should run its GPU path.

    False when device=cpu. When device=gpu: True if usable; otherwise raise
    :class:`GpuUnavailableError` (loud) unless ``PMA_DEVICE_FALLBACK=1`` is set, in
    which case warn + return False.
    """
    if requested_device() != "gpu":
        return False
    ok, reason = gpu_capable()
    if ok:
        return True
    if os.environ.get("PMA_DEVICE_FALLBACK", "0") == "1":
        _log(run_dir, stage,
             f"GPU requested but unavailable ({reason}); falling back to CPU "
             "because PMA_DEVICE_FALLBACK=1.", level="warning")
        return False
    raise GpuUnavailableError(
        f"compute.device=gpu requested but GPU is unavailable: {reason}. Fix the env/"
        "allocation (e.g. provision-env --device gpu, request --gres), or set "
        "PMA_DEVICE_FALLBACK=1 to explicitly allow a CPU fallback.")


def device_used(use: bool) -> str:
    return "gpu" if use else "cpu"


def _log(run_dir: str | Path | None, stage: str | None, message: str, *, level: str = "info") -> None:
    if run_dir is None:
        return
    try:
        from .log import log_event

        log_event(run_dir, {"stage": stage or "compute", "event": "device_dispatch",
                            "level": level, "message": message})
    except Exception:
        pass
