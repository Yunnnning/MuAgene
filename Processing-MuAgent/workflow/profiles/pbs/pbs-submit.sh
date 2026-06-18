#!/usr/bin/env bash
# pbs-submit.sh — invoked by snakemake's cluster-generic executor for each rule.
#
# Args (positional, in order; snakemake appends the jobscript path last):
#   $1  rule name
#   $2  snakemake-internal jobid
#   $3  threads (cpus_per_task)
#   $4  mem_mb (memory in megabytes)
#   $5  runtime in MINUTES (snakemake resource); converted to HH:MM:SS below
#   $6  jobscript path (appended by snakemake)
#
# Optional env vars:
#   PMA_PBS_QUEUE     → -q <queue>
#   PMA_PBS_PROJECT   → -P <project>
#   PMA_LOG_DIR       → log directory (default: ./logs)
#
# Preprocessing is CPU-only (_GPU_CAPABLE is empty in resources.smk). GPU routing
# for the integration subagent belongs in the integration pipeline's submit profile.
# PMA_DEVICE=cpu is always exported to child jobs, overriding any GPU value that
# the head-job may carry.
#
# PMA_SUBMIT_DRY_RUN=1 prints the resolved qsub command and exits 0 without
# submitting.
#
# Writes the PBS job id to stdout (qsub -terse), which snakemake parses to track
# the submitted job.
set -euo pipefail

if [ "$#" -lt 6 ]; then
    echo "pbs-submit.sh: expected 6 args (rule jobid threads mem_mb runtime jobscript), got $#" >&2
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

# Convert minutes → HH:MM:SS for PBS walltime.
h=$(( runtime_min / 60 ))
m=$(( runtime_min % 60 ))
walltime=$(printf '%02d:%02d:00' "$h" "$m")

log_dir="${PMA_LOG_DIR:-logs}"
mkdir -p "$log_dir"

# Preprocessing is always CPU. Override PMA_DEVICE=cpu so child jobs do not
# inherit PMA_DEVICE=gpu from a head job configured with --device gpu.
select_chunk="select=1:ncpus=${threads}:mem=${mem_mb}mb"
queue="${PMA_PBS_QUEUE:-}"
declare -a env_opts=(-v "PMA_DEVICE=cpu")

declare -a opts=(
    -terse
    -l "${select_chunk}"
    -l "walltime=${walltime}"
    -N "pma_${rule}"
    -j oe
    -o "${log_dir}/${rule}.${jobid}.log"
)
[ -n "${queue}" ] && opts+=(-q "$queue")
[ -n "${PMA_PBS_PROJECT:-}" ] && opts+=(-P "$PMA_PBS_PROJECT")
opts+=("${env_opts[@]}")

if [ "${PMA_SUBMIT_DRY_RUN:-0}" != "0" ]; then
    echo "DRY_RUN qsub ${opts[*]} ${jobscript}"
    exit 0
fi

# qsub -terse outputs just the job id (e.g. "1234567.pbs"); snakemake reads it as
# the submission handle.
exec qsub "${opts[@]}" "$jobscript"
