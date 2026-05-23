#!/usr/bin/env bash
#SBATCH --job-name=pma_runner
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=pma_runner-%j.out
#
# Head-job: runs the snakemake orchestrator that submits each per-stage child
# job to SLURM. Exits when the DAG completes or stops at a missing-approval gate.
#
# Submitted via `processing-muagent submit --executor slurm ...`, which sets the
# following env vars via `sbatch --export`:
#   PMA_CONFIG      path to the run.yaml configfile
#   PMA_TARGET      snakemake target (default: all)
#   PMA_REPO_ROOT   absolute path to processing-MuAgent/
# Plus any PMA_SLURM_* and PMA_NOTIFY_EMAIL vars from the submitter's env.
set -euo pipefail

cd "${PMA_REPO_ROOT}"

SNAKE_ARGS=(
    --profile "workflow/profiles/slurm"
    --configfile "${PMA_CONFIG}"
    "${PMA_TARGET:-all}"
)

EXIT_CODE=0
scripts/launch_runner.sh "${SNAKE_ARGS[@]}" || EXIT_CODE=$?

if [ -n "${PMA_NOTIFY_EMAIL:-}" ]; then
    subject="processing-MuAgent: runner exited with code ${EXIT_CODE}"
    body="Head-job ${SLURM_JOB_ID:-unknown} finished.
Config:   ${PMA_CONFIG}
Target:   ${PMA_TARGET:-all}
ExitCode: ${EXIT_CODE}

If exit code is 0 and target was 'all', the run is complete.
If exit code is 0 and the runner stopped early, snakemake hit an
unapproved checkpoint (typically s7_clustering). Review the proposal
artifact at <run_dir>/internal/proposals/<stage>.yaml, then
'processing-muagent approve <stage>' and 'processing-muagent submit'
again to resume.

Run dir: $(awk '/^run_dir:/ {print $2}' "${PMA_CONFIG}")
"
    if command -v mail >/dev/null 2>&1; then
        printf '%s\n' "$body" | mail -s "$subject" "$PMA_NOTIFY_EMAIL" || true
    elif command -v sendmail >/dev/null 2>&1; then
        { printf 'To: %s\nSubject: %s\n\n%s\n' "$PMA_NOTIFY_EMAIL" "$subject" "$body"; } \
            | sendmail -t || true
    fi
fi

exit "${EXIT_CODE}"
