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

# Put the repo on PYTHONPATH so `import executor` resolves from source in EVERY job tier —
# this process AND the child jobs snakemake submits (they inherit it via sbatch/qsub
# --export=ALL). The head job already finds executor via cwd=repo_root and GPU child jobs
# via the container wrapper's --env PYTHONPATH; this closes the CPU child-job gap so a
# submit-time auto-provisioned env (created from the lock, no `pip install -e`) still
# imports executor. init-machine's editable install remains for interactive console scripts.
export PYTHONPATH="${PMA_REPO_ROOT}:${PYTHONPATH:-}"

# Activate the project conda env. Identity comes from PMA_CONDA_ENV (set by
# configure-execution); `muagene` is the canonical default for a fresh install.
CONDA_ENV="${PMA_CONDA_ENV:-muagene}"
if [ -z "${PMA_CONDA_ENV:-}" ]; then
    echo "launch_runner.sh: PMA_CONDA_ENV unset; defaulting to '$CONDA_ENV'. Set it with" \
         "configure-execution --conda-env <name> (or provision it via Execution-MuAgent provision-env)." >&2
fi

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
if [ -n "$configfile" ]; then
    run_dir="$(python - "$configfile" <<'PY'
import sys
from pathlib import Path
import yaml

with Path(sys.argv[1]).open() as fh:
    cfg = yaml.safe_load(fh) or {}
print(Path(cfg["run_dir"]).expanduser().resolve())
PY
)"
    # Single live source for the run directory: export it so the child-job submit
    # scripts (spawned by snakemake under this process) can bind it into GPU
    # containers. Derived from the config here — deliberately NOT projected into
    # hpc.env, which is a static site.config snapshot a copied path would drift from.
    export PMA_RUN_DIR="$run_dir"
    if [ "$has_directory_arg" -eq 0 ]; then
        snakemake_workdir="${run_dir}/internal/snakemake"
        mkdir -p "$snakemake_workdir"
        export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${snakemake_workdir}/cache}"
        directory_args=(--directory "$snakemake_workdir")
    fi
fi

# Cluster profiles: shared-NFS snakemake flags + site-specific SLURM defaults.
extra_args=()
using_cluster_profile=0
for arg in "$@"; do
    if [[ "$arg" == */profiles/slurm* ]] || [[ "$arg" == */profiles/pbs* ]]; then
        using_cluster_profile=1
    fi
    if [[ "$arg" == */profiles/slurm* ]]; then
        if [ -n "${PMA_SLURM_PARTITION:-}" ]; then
            extra_args+=(--default-resources "slurm_partition=$PMA_SLURM_PARTITION")
        fi
        if [ -n "${PMA_SLURM_ACCOUNT:-}" ]; then
            extra_args+=(--default-resources "slurm_account=$PMA_SLURM_ACCOUNT")
        fi
    fi
done
if [ "$using_cluster_profile" -eq 1 ]; then
    extra_args+=(
        --shared-fs-usage persistence input-output software-deployment
        software-deployment-cache sources source-cache
    )
fi

exec python -m snakemake -s "$REPO_ROOT/workflow/Snakefile" "${directory_args[@]}" "$@" "${extra_args[@]}"
