#!/usr/bin/env bash
# slurm-status.sh - invoked by Snakemake's cluster-generic executor to poll jobs.
#
# Arg: <slurm_job_id>
# Output (one of): success | failed | running
set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "slurm-status.sh: expected 1 arg (slurm_job_id), got $#" >&2
    exit 2
fi

jobid="${1%%;*}"
query_timeout="${PMA_SCHEDULER_QUERY_TIMEOUT:-5}"

state=""
if command -v sacct >/dev/null 2>&1; then
    state=$(timeout "$query_timeout" sacct -j "$jobid" -X -n -P -o State 2>/dev/null \
        | awk -F'|' 'NF {print $1; exit}') || true
fi

if [ -z "$state" ] && command -v squeue >/dev/null 2>&1; then
    state=$(timeout "$query_timeout" squeue -j "$jobid" -h -o "%T" 2>/dev/null | awk 'NF {print $1; exit}') || true
fi

case "$state" in
    COMPLETED)
        echo success
        ;;
    PENDING|RUNNING|CONFIGURING|COMPLETING|SUSPENDED|RESIZING)
        echo running
        ;;
    "")
        # Scheduler accounting can lag or hang on busy clusters. Unknown state
        # is treated as running so Snakemake does not falsely accept a failed job
        # before declared outputs exist.
        echo running
        ;;
    *)
        echo failed
        ;;
esac
