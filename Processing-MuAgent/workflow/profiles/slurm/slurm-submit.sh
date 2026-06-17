#!/usr/bin/env bash
# slurm-submit.sh - invoked by Snakemake's cluster-generic executor.
#
# Args (positional, in order; Snakemake appends the jobscript path last):
#   $1  rule name
#   $2  snakemake-internal jobid
#   $3  threads (cpus-per-task)
#   $4  mem_mb (memory in megabytes)
#   $5  runtime in minutes
#   $6  gpu count (0 for CPU rules; >0 for GPU-capable rules — see resources.smk)
#   $7  jobscript path
#
# GPU routing (when $6 > 0): land on $PMA_SLURM_GPU_PARTITION (falls back to the
# normal partition) with --gres=$PMA_SLURM_GPU_GRES, and provide the GPU env per
# provider:
#   PMA_GPU_PROVIDER=container -> run the jobscript inside `singularity exec --nv
#                                 $PMA_GPU_IMAGE` (the image carries python+rapids)
#   else (conda-env provider)  -> export PMA_CONDA_ENV=$PMA_CONDA_ENV_GPU so the
#                                 injected `conda activate` brings up the GPU env
#
# PMA_SUBMIT_DRY_RUN=1 prints the resolved sbatch command (and any container
# wrapper) and exits 0 without submitting — used by render/parity tests.
set -euo pipefail

if [ "$#" -lt 7 ]; then
    echo "slurm-submit.sh: expected 7 args (rule jobid threads mem_mb runtime gpu jobscript), got $#" >&2
    exit 2
fi

rule="$1"
jobid="$2"
threads="$3"
mem_mb="$4"
runtime_min="$5"
gpu="${6:-0}"
jobscript="$7"

if [ "${PMA_DISABLE_STORAGE_LOCAL_COPIES:-1}" != "0" ]; then
    PYTHONPATH="${PMA_REPO_ROOT}:${PYTHONPATH:-}" python - "$jobscript" <<'PY'
import sys
from executor import hpc

hpc.sanitize_snakemake_jobscript(sys.argv[1])
PY
fi

log_root="${PMA_LOG_DIR:-.snakemake/slurm_logs}"
log_dir="${log_root}/rule_${rule}"
mkdir -p "$log_dir"

# Defaults preserve the existing CPU behaviour exactly: --export=ALL on the normal
# partition. GPU-capable jobs override partition/gres/env below.
export_spec="ALL"
partition="${PMA_SLURM_PARTITION:-}"
declare -a gpu_opts=()
is_gpu=0
if [ "${gpu}" != "0" ] && [ -n "${gpu}" ]; then
    is_gpu=1
    [ -n "${PMA_SLURM_GPU_PARTITION:-}" ] && partition="$PMA_SLURM_GPU_PARTITION"
    [ -n "${PMA_SLURM_GPU_GRES:-}" ] && gpu_opts+=(--gres "$PMA_SLURM_GPU_GRES")
    if [ "${PMA_GPU_PROVIDER:-}" = "container" ]; then
        export_spec="ALL,PMA_DEVICE=gpu"
    else
        # conda-env GPU provider: override PMA_CONDA_ENV so the injected activation
        # brings up the GPU env instead of the inherited CPU env.
        export_spec="ALL,PMA_DEVICE=gpu${PMA_CONDA_ENV_GPU:+,PMA_CONDA_ENV=${PMA_CONDA_ENV_GPU}}"
    fi
fi

declare -a opts=(
    --parsable
    --export="${export_spec}"
    --cpus-per-task="${threads}"
    --mem="${mem_mb}M"
    --time="${runtime_min}"
    --job-name="pma_${rule}"
    --output="${log_dir}/%j.log"
)
[ -n "${partition}" ] && opts+=(--partition "$partition")
[ -n "${PMA_SLURM_ACCOUNT:-}" ] && opts+=(--account "$PMA_SLURM_ACCOUNT")
[ "${#gpu_opts[@]}" -gt 0 ] && opts+=("${gpu_opts[@]}")

# GPU container bind contract — KEEP IDENTICAL in pbs-submit.sh:
#   A container that runs pipeline code MUST bind BOTH the resolved repo root
#   (PMA_REPO_ROOT — for launch_runner.sh + the `executor` package on PYTHONPATH)
#   AND the resolved run directory (PMA_RUN_DIR — for internal/artifacts I/O),
#   because the run data may sit under a nested mount that singularity's default
#   $HOME/$PWD auto-mount does not cover. Optional extra binds (PMA_GPU_BIND,
#   sourced from site.config common.scratch via hpc.env) append after. Producers
#   resolve these paths (launch_runner.sh / configure-execution), so no re-resolve here.
# For a GPU container job, submit a thin wrapper that exec's the jobscript inside
# the image with GPU passthrough; otherwise submit the jobscript directly.
target="$jobscript"
if [ "${is_gpu}" = "1" ] && [ "${PMA_GPU_PROVIDER:-}" = "container" ] && [ -n "${PMA_GPU_IMAGE:-}" ]; then
    if [ -z "${PMA_RUN_DIR:-}" ]; then
        echo "slurm-submit.sh: WARNING — GPU container job but PMA_RUN_DIR is unset; the run" \
             "directory will not be bound and the job cannot read inputs / write outputs." \
             "Ensure launch_runner.sh exported it (a configfile must be on the snakemake CLI)." >&2
    fi
    target="${jobscript}.gpuwrap.sh"
    {
        echo "#!/usr/bin/env bash"
        echo "set -euo pipefail"
        [ -n "${PMA_SINGULARITY_MODULE:-}" ] && echo "module load ${PMA_SINGULARITY_MODULE} 2>/dev/null || true"
        printf 'exec singularity exec --nv'
        printf ' --env PYTHONPATH=%q' "${PMA_REPO_ROOT}"
        printf ' --bind %q' "${PMA_REPO_ROOT}"
        [ -n "${PMA_RUN_DIR:-}" ] && printf ' --bind %q' "${PMA_RUN_DIR}"
        [ -n "${PMA_GPU_BIND:-}" ] && printf ' --bind %q' "${PMA_GPU_BIND}"
        printf ' %q bash %q\n' "${PMA_GPU_IMAGE}" "${jobscript}"
    } > "$target"
    chmod +x "$target"
fi

if [ "${PMA_SUBMIT_DRY_RUN:-0}" != "0" ]; then
    echo "DRY_RUN sbatch ${opts[*]} ${target}"
    if [ "$target" != "$jobscript" ]; then
        echo "--- container wrapper ($target) ---"
        cat "$target"
    fi
    exit 0
fi

exec sbatch "${opts[@]}" "$target"
