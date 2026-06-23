"""Per-stage resource declarations — single source of truth for SLURM profiles.

Included by the top-level Snakefile. The values in `RESOURCES` and `RUNTIME` are
referenced by `<stage>_execute` rules via `resources: mem_mb=..., runtime=..., cpus_per_task=...`.

Conventions:
- `mem_mb` is total job memory in megabytes (snakemake-standard resource name).
- `runtime` is wall-clock walltime in MINUTES (snakemake-standard; required by
  snakemake>=8). SLURM accepts minutes directly.
- `cpus` becomes `cpus_per_task` (SLURM) via the profile templates.

Defaults are sized for ~10 000-cell datasets. Two scaling knobs:
  1. PMA_RESOURCES_SCALE environment variable: multiplies mem_mb and walltime
     (e.g. PMA_RESOURCES_SCALE=2 for ~30 k-cell, =4 for ~100 k-cell).
  2. Snakemake `restart-times` + `attempt`-aware memory in each rule
     (defined in the .smk files): a job killed for OOM is resubmitted with 2x mem.
"""
from __future__ import annotations

import os


def _scale_factor() -> float:
    try:
        return max(1.0, float(os.environ.get("PMA_RESOURCES_SCALE", "1")))
    except ValueError:
        return 1.0


def _scaled_mem(mem_mb: int) -> int:
    return int(mem_mb * _scale_factor())


def _scaled_runtime_min(minutes: int) -> int:
    """Scale a walltime (minutes) by PMA_RESOURCES_SCALE and add an NFS overhead buffer.

    Shared network filesystems can add slow post-job visibility/write-back
    overhead around large h5ad outputs. Snakemake 9 cluster child jobs also
    hang if storage-local-copies is enabled on NFS (see workflow/profiles/). The buffer absorbs this cost without
    setting PMA_RESOURCES_SCALE. PMA_RESOURCES_SCALE still multiplies the base
    compute time; the buffer is added afterwards so large-dataset scaling
    doesn't double the buffer.
    """
    _NFS_OVERHEAD_MIN = int(os.environ.get("PMA_RUNTIME_OVERHEAD_MIN", "90"))
    return int(minutes * _scale_factor() + _NFS_OVERHEAD_MIN)


# Base resource table. Edit here; the SLURM profile picks this up.
_BASE_RESOURCES: dict[str, dict[str, int]] = {
    # Local rules — declared for completeness.
    "p1_context":     {"cpus": 1, "mem_mb": 1_000},
    "plan_review":    {"cpus": 1, "mem_mb": 1_000},
    # s0_ingest is the merged planning compute: load RNA + import ATAC fragments +
    # full RNA QC + threshold exploration in one cluster job — sized like s2_atac_qc.
    "s0_ingest":      {"cpus": 4, "mem_mb": 32_000},
    # Cluster rules.
    "s1a_ambient":    {"cpus": 2, "mem_mb": 16_000},
    "s1_rna_qc":      {"cpus": 1, "mem_mb": 8_000},
    "s2_atac_qc":     {"cpus": 4, "mem_mb": 32_000},
    "s3_doublets":    {"cpus": 2, "mem_mb": 32_000},
    "s4_rna_norm":    {"cpus": 1, "mem_mb": 8_000},
    "s5_atac_spectral":    {"cpus": 4, "mem_mb": 32_000},
    "s6_neighbors":      {"cpus": 2, "mem_mb": 16_000},
    "s7_clustering":  {"cpus": 4, "mem_mb": 16_000},
    "s8_umap":        {"cpus": 2, "mem_mb": 8_000},
    # s_handoff is a localrule (orchestrator host): assemble post-QC .h5mu + manifest.
    "s_handoff":      {"cpus": 1, "mem_mb": 4_000},
}

# Walltime in MINUTES. snakemake>=8 requires `runtime` to be int minutes.
_BASE_RUNTIME_MIN: dict[str, int] = {
    "p1_context":     15,
    "plan_review":    10,
    "s0_ingest":     360,
    "s1a_ambient":    60,
    "s1_rna_qc":     120,
    "s2_atac_qc":    360,
    "s3_doublets":   240,
    "s4_rna_norm":    30,
    "s5_atac_spectral":    90,
    "s6_neighbors":      45,
    "s7_clustering": 120,
    "s8_umap":        45,
    "s_handoff":      30,
}


# GPU routing ----------------------------------------------------------------
# _GPU_CAPABLE is the registry of stages that can use a GPU when compute.device=gpu.
# Preprocessing is CPU-only — the set is empty. GPU routing (partition, gres,
# container dispatch) lives in the integration pipeline's submit profile, not here.
#
# When a preprocessing stage gains GPU support in the future, three edits are
# needed in lockstep:
#   1. Add the stage name to _GPU_CAPABLE (the gate).
#   2. Re-wire the submit profile (slurm/config.yaml + slurm-submit.sh) to
#      accept and route a `gpu` resource — the preprocessing profiles intentionally
#      omit {resources.gpu} today because no preprocessing stage uses it.
#   3. Give the stage a `compute.use_gpu()` branch for its compute path.
_GPU_CAPABLE: set[str] = set()


def _device() -> str:
    return (os.environ.get("PMA_DEVICE", "cpu") or "cpu").strip().lower()


def gpu_for(stage: str) -> int:
    """GPUs to request for a stage: 1 when device=gpu and the stage is GPU-capable, else 0."""
    return 1 if (_device() == "gpu" and stage in _GPU_CAPABLE) else 0


# Public API ----------------------------------------------------------------

# `gpu` is intentionally absent from RESOURCES: preprocessing submit profiles do not
# pass a gpu resource to the scheduler (see slurm/config.yaml). When a stage is
# added to _GPU_CAPABLE in the future, the submit profiles must be re-wired first.
RESOURCES: dict[str, dict[str, int]] = {
    name: {"cpus": v["cpus"], "mem_mb": _scaled_mem(v["mem_mb"])}
    for name, v in _BASE_RESOURCES.items()
}

# Walltime in MINUTES — directly usable as a snakemake `runtime` resource.
RUNTIME: dict[str, int] = {name: _scaled_runtime_min(v) for name, v in _BASE_RUNTIME_MIN.items()}


def mem_mb_for(stage: str, attempt: int = 1) -> int:
    """Return mem_mb for a given attempt (1-based). Doubles on retry to absorb OOM kills.

    Use in rule resources as:
        resources: mem_mb=lambda wc, attempt: mem_mb_for("s3_doublets", attempt)
    """
    return RESOURCES[stage]["mem_mb"] * attempt


# Progress-timeout hints in MINUTES — consumed by Execution-MuAgent monitor.
# Single source of truth for per-stage silence thresholds. specs.py reads this dict
# and writes the values into internal/stage_meta/<stage>.yaml at plan-review time.
# Execution-MuAgent reads those YAMLs; --stale-minutes 90 is the fallback default
# when no hint is present (e.g. for the head_job itself).
# Independent of PMA_RESOURCES_SCALE: reflects algorithm cadence, not compute time.
PROGRESS_TIMEOUT_HINT: dict[str, int] = {
    "s0_ingest":       120,
    "s1a_ambient":      30,
    "s1_rna_qc":        45,
    "s2_atac_qc":      120,
    "s3_doublets":      60,
    "s4_rna_norm":      20,
    "s5_atac_spectral": 45,
    "s6_neighbors":     30,
    "s7_clustering":    60,
    "s8_umap":          30,
    "s_handoff":        20,
}
