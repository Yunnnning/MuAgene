#!/usr/bin/env bash
# pbs-status.sh — invoked by snakemake's cluster-generic executor to poll job state.
#
# Arg: <pbs_job_id>
# Output (one of): success | failed | running
set -euo pipefail

if [ "$#" -lt 1 ]; then
    echo "pbs-status.sh: expected 1 arg (pbs_job_id), got $#" >&2
    exit 2
fi

jobid="$1"
query_timeout="${PMA_SCHEDULER_QUERY_TIMEOUT:-5}"

# Try the live queue first; if absent, try the finished-job history (-x).
state=""
state=$(timeout "$query_timeout" qstat -f "$jobid" 2>/dev/null | awk -F= '/job_state/ {gsub(/ /, "", $2); print $2; exit}') || true
if [ -z "$state" ]; then
    state=$(timeout "$query_timeout" qstat -fx "$jobid" 2>/dev/null | awk -F= '/job_state/ {gsub(/ /, "", $2); print $2; exit}') || true
fi

case "$state" in
    R|Q|H|S|T|W|E|B|M)
        echo running
        ;;
    F)
        # Finished — check Exit_status (only present after job completes).
        ec=$(timeout "$query_timeout" qstat -fx "$jobid" 2>/dev/null | awk -F= '/Exit_status/ {gsub(/ /, "", $2); print $2; exit}') || true
        if [ "${ec:-0}" = "0" ]; then
            echo success
        else
            echo failed
        fi
        ;;
    "")
        # Unknown state can mean scheduler lag or a timed-out query. Keep
        # polling rather than falsely accepting a failed or cancelled job.
        echo running
        ;;
    *)
        echo failed
        ;;
esac
