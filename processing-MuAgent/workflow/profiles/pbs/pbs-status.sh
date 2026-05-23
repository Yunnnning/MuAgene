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

# Try the live queue first; if absent, try the finished-job history (-x).
state=""
state=$(qstat -f "$jobid" 2>/dev/null | awk -F= '/job_state/ {gsub(/ /, "", $2); print $2; exit}') || true
if [ -z "$state" ]; then
    state=$(qstat -fx "$jobid" 2>/dev/null | awk -F= '/job_state/ {gsub(/ /, "", $2); print $2; exit}') || true
fi

case "$state" in
    R|Q|H|S|T|W|E|B|M)
        echo running
        ;;
    F)
        # Finished — check Exit_status (only present after job completes).
        ec=$(qstat -fx "$jobid" 2>/dev/null | awk -F= '/Exit_status/ {gsub(/ /, "", $2); print $2; exit}') || true
        if [ "${ec:-0}" = "0" ]; then
            echo success
        else
            echo failed
        fi
        ;;
    "")
        # Job not visible in either live queue or history → assume it finished
        # cleanly and was purged from history (some sites have short retention).
        # Returning "failed" here would cause snakemake to mark the rule failed
        # even though the output may exist; let snakemake decide by checking the
        # output files. "success" is the safer default when state is unknown.
        echo success
        ;;
    *)
        echo failed
        ;;
esac
