#!/usr/bin/env bash
# slurm-submit.sh - invoked by Snakemake's cluster-generic executor.
#
# Args (positional, in order; Snakemake appends the jobscript path):
#   $1  rule name
#   $2  snakemake-internal jobid
#   $3  threads (cpus-per-task)
#   $4  mem_mb (memory in megabytes)
#   $5  runtime in minutes
#   $6  jobscript path
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

log_root="${PMA_LOG_DIR:-.snakemake/slurm_logs}"
log_dir="${log_root}/rule_${rule}"
mkdir -p "$log_dir"

declare -a opts=(
    --parsable
    --export=ALL
    --cpus-per-task="${threads}"
    --mem="${mem_mb}M"
    --time="${runtime_min}"
    --job-name="pma_${rule}"
    --output="${log_dir}/%j.log"
)

if [ -n "${PMA_SLURM_PARTITION:-}" ]; then
    opts+=(--partition "$PMA_SLURM_PARTITION")
fi

if [ -n "${PMA_SLURM_ACCOUNT:-}" ]; then
    opts+=(--account "$PMA_SLURM_ACCOUNT")
fi

exec sbatch "${opts[@]}" "$jobscript"
