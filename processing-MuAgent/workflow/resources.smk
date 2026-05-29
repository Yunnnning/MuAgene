"""Per-stage resource declarations — single source of truth for PBS and SLURM profiles.

Included by the top-level Snakefile. The values in `RESOURCES` and `RUNTIME` are
referenced by `<stage>_execute` rules via `resources: mem_mb=..., runtime=..., cpus_per_task=...`.

Conventions:
- `mem_mb` is total job memory in megabytes (snakemake-standard resource name).
- `runtime` is wall-clock walltime in MINUTES (snakemake-standard; required by
  snakemake>=8). PBS profile converts to `HH:MM:SS` in its submit script;
  SLURM accepts minutes directly.
- `cpus` becomes `cpus_per_task` (SLURM) / `ncpus` (PBS) via the profile templates.

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

    On shared network filesystems Snakemake 9's post-job storage scan can add
    ~30 min of overhead per child job even when no files need to be transferred
    (snakemake-executor-plugin-slurm-jobstep + store_storage_outputs on NFS).
    A 60-min buffer absorbs this cost without setting PMA_RESOURCES_SCALE.
    PMA_RESOURCES_SCALE still multiplies the *base* compute time; the buffer
    is added afterwards so large-dataset scaling doesn't double the buffer.
    """
    _NFS_OVERHEAD_MIN = int(os.environ.get("PMA_RUNTIME_OVERHEAD_MIN", "60"))
    return int(minutes * _scale_factor() + _NFS_OVERHEAD_MIN)


# Base resource table. Edit here; both PBS and SLURM profiles pick this up.
_BASE_RESOURCES: dict[str, dict[str, int]] = {
    # Local rules — declared for completeness; s0_ingest_execute uses these on cluster.
    "p1_context":     {"cpus": 1, "mem_mb": 1_000},
    "p2_plan":        {"cpus": 1, "mem_mb": 1_000},
    "plan_review":    {"cpus": 1, "mem_mb": 1_000},
    "s0_ingest":      {"cpus": 1, "mem_mb": 8_000},
    # Cluster rules.
    "s1a_ambient":    {"cpus": 2, "mem_mb": 16_000},
    "s1_rna_qc":      {"cpus": 1, "mem_mb": 8_000},
    "s2_atac_qc":     {"cpus": 2, "mem_mb": 16_000},
    "s3_doublets":    {"cpus": 2, "mem_mb": 32_000},
    "s4_rna_norm":    {"cpus": 1, "mem_mb": 8_000},
    "s5_atac_lsi":    {"cpus": 4, "mem_mb": 32_000},
    "s6_dimred":      {"cpus": 2, "mem_mb": 16_000},
    "s7_clustering":  {"cpus": 4, "mem_mb": 16_000},
    "s8_umap":        {"cpus": 2, "mem_mb": 8_000},
}

# Walltime in MINUTES. snakemake>=8 requires `runtime` to be int minutes.
_BASE_RUNTIME_MIN: dict[str, int] = {
    "p1_context":     15,
    "p2_plan":        15,
    "plan_review":    10,
    "s0_ingest":      30,
    "s1a_ambient":    60,
    "s1_rna_qc":      30,
    "s2_atac_qc":     60,
    "s3_doublets":   120,
    "s4_rna_norm":    30,
    "s5_atac_lsi":   120,
    "s6_dimred":      45,
    "s7_clustering": 120,
    "s8_umap":        45,
}


# Public API ----------------------------------------------------------------

RESOURCES: dict[str, dict[str, int]] = {
    name: {"cpus": v["cpus"], "mem_mb": _scaled_mem(v["mem_mb"])}
    for name, v in _BASE_RESOURCES.items()
}

# Walltime in MINUTES — directly usable as a snakemake `runtime` resource.
# PBS submit script converts minutes → HH:MM:SS at qsub time.
RUNTIME: dict[str, int] = {name: _scaled_runtime_min(v) for name, v in _BASE_RUNTIME_MIN.items()}


def mem_mb_for(stage: str, attempt: int = 1) -> int:
    """Return mem_mb for a given attempt (1-based). Doubles on retry to absorb OOM kills.

    Use in rule resources as:
        resources: mem_mb=lambda wc, attempt: mem_mb_for("s3_doublets", attempt)
    """
    return RESOURCES[stage]["mem_mb"] * attempt
