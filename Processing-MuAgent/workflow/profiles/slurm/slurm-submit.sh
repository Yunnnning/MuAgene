#!/usr/bin/env bash
# slurm-submit.sh - invoked by Snakemake's cluster-generic executor.
#
# Args (positional, in order; Snakemake appends the jobscript path last):
#   $1  rule name
#   $2  snakemake-internal jobid
#   $3  threads (cpus-per-task)
#   $4  mem_mb (memory in megabytes)
#   $5  runtime in minutes
#   $6  jobscript path
#
# Preprocessing is CPU-only (_GPU_CAPABLE is empty in resources.smk). GPU routing
# for the integration subagent belongs in the integration pipeline's submit profile.
# PMA_DEVICE=cpu is always exported to child jobs, overriding any GPU value that
# the head-job may carry (head job is configured with --device gpu for integration,
# but preprocessing child jobs must never inherit it).
#
# PMA_SUBMIT_DRY_RUN=1 prints the resolved sbatch command and exits 0 without
# submitting — used by render/parity tests.
set -euo pipefail

if [ "$#" -lt 6 ]; then
    echo "slurm-submit.sh: expected 6 args (rule jobid threads mem_mb runtime jobscript), got $#" >&2
    exit 2
fi

rule="$1"
jobid="$2"
threads="$3"
mem_mb="$4"
runtime_min="$5"
jobscript="$6"

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

# Preprocessing is always CPU. Override PMA_DEVICE=cpu so child jobs do not
# inherit PMA_DEVICE=gpu from a head job configured with --device gpu.
export_spec="ALL,PMA_DEVICE=cpu"
partition="${PMA_SLURM_PARTITION:-}"

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

if [ "${PMA_SUBMIT_DRY_RUN:-0}" != "0" ]; then
    echo "DRY_RUN sbatch ${opts[*]} ${jobscript}"
    exit 0
fi

exec sbatch "${opts[@]}" "$jobscript"
