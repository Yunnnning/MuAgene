#!/usr/bin/env bash
#SBATCH --job-name=pma_runner
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=logs/pma_runner-%j.out
#
# Head-job: runs the snakemake orchestrator that submits each per-stage child
# job to SLURM. Exits when the DAG completes or stops at a missing-approval gate.
#
# Submitted by Execution-MuAgent (execute-spec), which renders this script from
# the head_job.yaml spec and sets the following env vars:
#   PMA_CONFIG      path to the run.yaml configfile
#   PMA_TARGET      snakemake target (default: all)
#   PMA_REPO_ROOT   absolute path to Processing-MuAgent/
# Monitoring and findings are written to internal/hpc_monitor/ by Execution-MuAgent.
set -euo pipefail

cd "${PMA_REPO_ROOT}"

SNAKE_ARGS=(
    --profile "workflow/profiles/slurm"
    "${PMA_TARGET:-all}"
    --configfile "${PMA_CONFIG}"
)

scripts/launch_runner.sh "${SNAKE_ARGS[@]}"
