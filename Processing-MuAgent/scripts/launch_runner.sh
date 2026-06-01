#!/usr/bin/env bash
# launch_runner.sh — activate project env and invoke snakemake.
#
# Used by both interactive mode (called from executor.cli._snakemake when
# --executor != local) and headless mode (called from runner.pbs / runner.slurm).
#
# All trailing args are forwarded to snakemake. Typical invocation:
#   launch_runner.sh --profile workflow/profiles/pbs --configfile $CFG all
set -euo pipefail

# Locate the Processing-MuAgent repo root (parent of this script's dir).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"
export PMA_REPO_ROOT="${PMA_REPO_ROOT:-$REPO_ROOT}"

# Activate the project conda env. Allow override via PMA_CONDA_ENV.
CONDA_ENV="${PMA_CONDA_ENV:-grn}"

if command -v micromamba >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate "$CONDA_ENV"
elif command -v mamba >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    eval "$(mamba shell hook --shell bash 2>/dev/null || conda shell.bash hook)"
    mamba activate "$CONDA_ENV" || conda activate "$CONDA_ENV"
elif command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV"
else
    echo "launch_runner.sh: no conda/mamba/micromamba on PATH; cannot activate $CONDA_ENV" >&2
    exit 2
fi

# Preserve the reproducibility-oriented thread caps the CLI sets in local mode.
# On HPC the user can opt into more threads by setting OMP_NUM_THREADS in their
# job script, but the default mirrors local behaviour.
export NUMBA_NUM_THREADS="${NUMBA_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export PYTHONHASHSEED="${PYTHONHASHSEED:-0}"

configfile=""
has_directory_arg=0
prev=""
for arg in "$@"; do
    if [ "$prev" = "--configfile" ]; then
        configfile="$arg"
    fi
    if [ "$arg" = "--directory" ] || [[ "$arg" == --directory=* ]]; then
        has_directory_arg=1
    fi
    prev="$arg"
done

directory_args=()
if [ "$has_directory_arg" -eq 0 ] && [ -n "$configfile" ]; then
    run_dir="$(python - "$configfile" <<'PY'
import sys
from pathlib import Path
import yaml

with Path(sys.argv[1]).open() as fh:
    cfg = yaml.safe_load(fh) or {}
print(Path(cfg["run_dir"]).expanduser().resolve())
PY
)"
    snakemake_workdir="${run_dir}/internal/snakemake"
    mkdir -p "$snakemake_workdir"
    export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${snakemake_workdir}/cache}"
    directory_args=(--directory "$snakemake_workdir")
fi

# Add site-specific --default-resources for SLURM partition / account if the
# user has set the env vars. Detected from the --profile arg below.
extra_args=()
for arg in "$@"; do
    if [[ "$arg" == */profiles/slurm* ]]; then
        if [ -n "${PMA_SLURM_PARTITION:-}" ]; then
            extra_args+=(--default-resources "slurm_partition=$PMA_SLURM_PARTITION")
        fi
        if [ -n "${PMA_SLURM_ACCOUNT:-}" ]; then
            extra_args+=(--default-resources "slurm_account=$PMA_SLURM_ACCOUNT")
        fi
        break
    fi
done

exec python -m snakemake -s "$REPO_ROOT/workflow/Snakefile" "${directory_args[@]}" "$@" "${extra_args[@]}"
