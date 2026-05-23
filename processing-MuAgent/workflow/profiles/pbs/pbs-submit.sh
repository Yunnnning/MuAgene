#!/usr/bin/env bash
# pbs-submit.sh — invoked by snakemake's cluster-generic executor for each rule.
#
# Args (positional, in order; snakemake appends the jobscript path):
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
# Writes the PBS job id to stdout (qsub -terse), which snakemake parses to track
# the submitted job.
set -euo pipefail

if [ "$#" -lt 6 ]; then
    echo "pbs-submit.sh: expected 6 args (rule jobid threads mem_mb walltime jobscript), got $#" >&2
    exit 2
fi

rule="$1"
jobid="$2"
threads="$3"
mem_mb="$4"
runtime_min="$5"
jobscript="$6"

# Convert minutes → HH:MM:SS for PBS walltime.
h=$(( runtime_min / 60 ))
m=$(( runtime_min % 60 ))
walltime=$(printf '%02d:%02d:00' "$h" "$m")

log_dir="${PMA_LOG_DIR:-logs}"
mkdir -p "$log_dir"

declare -a opts=(
    -terse
    -l "select=1:ncpus=${threads}:mem=${mem_mb}mb"
    -l "walltime=${walltime}"
    -N "pma_${rule}"
    -j oe
    -o "${log_dir}/${rule}.${jobid}.log"
)

if [ -n "${PMA_PBS_QUEUE:-}" ]; then
    opts+=(-q "$PMA_PBS_QUEUE")
fi

if [ -n "${PMA_PBS_PROJECT:-}" ]; then
    opts+=(-P "$PMA_PBS_PROJECT")
fi

# qsub -terse outputs just the job id (e.g. "1234567.pbs"); snakemake reads it as
# the submission handle.
exec qsub "${opts[@]}" "$jobscript"
