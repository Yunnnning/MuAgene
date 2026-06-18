#!/usr/bin/env bash
# pbs-submit.sh — invoked by snakemake's cluster-generic executor for each rule.
#
# Args (positional, in order; snakemake appends the jobscript path last):
#   $1  rule name
#   $2  snakemake-internal jobid
#   $3  threads (cpus_per_task)
#   $4  mem_mb (memory in megabytes)
#   $5  runtime in MINUTES (snakemake resource); converted to HH:MM:SS below
#   $6  gpu count (0 for CPU rules; >0 for GPU-capable rules — see resources.smk)
#   $7  jobscript path (appended by snakemake)
#
# Optional env vars:
#   PMA_PBS_QUEUE     → -q <queue>
#   PMA_PBS_PROJECT   → -P <project>
#   PMA_LOG_DIR       → log directory (default: ./logs)
#
# GPU routing (when $6 > 0): append $PMA_PBS_GPU_SELECT_EXTRA (e.g. "ngpus=1" or
# "ngpus=1:gpu_type=a100" — PBS GPU syntax is site-variable, hence a template) to
# the select chunk; optionally route to $PMA_PBS_GPU_QUEUE; provide the GPU env per
# provider (container -> singularity exec --nv; else conda env via -v PMA_CONDA_ENV).
#
# PMA_SUBMIT_DRY_RUN=1 prints the resolved qsub command (and any container wrapper)
# and exits 0 without submitting.
#
# Writes the PBS job id to stdout (qsub -terse), which snakemake parses to track
# the submitted job.
set -euo pipefail

if [ "$#" -lt 7 ]; then
    echo "pbs-submit.sh: expected 7 args (rule jobid threads mem_mb runtime gpu jobscript), got $#" >&2
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

# Convert minutes → HH:MM:SS for PBS walltime.
h=$(( runtime_min / 60 ))
m=$(( runtime_min % 60 ))
walltime=$(printf '%02d:%02d:00' "$h" "$m")

log_dir="${PMA_LOG_DIR:-logs}"
mkdir -p "$log_dir"

# Base select chunk; GPU-capable jobs append the site's GPU template and may route
# to a dedicated GPU queue. Defaults preserve the existing CPU behaviour exactly.
select_chunk="select=1:ncpus=${threads}:mem=${mem_mb}mb"
queue="${PMA_PBS_QUEUE:-}"
declare -a env_opts=()
is_gpu=0
if [ "${gpu}" != "0" ] && [ -n "${gpu}" ]; then
    is_gpu=1
    [ -n "${PMA_PBS_GPU_SELECT_EXTRA:-}" ] && select_chunk="${select_chunk}:${PMA_PBS_GPU_SELECT_EXTRA}"
    [ -n "${PMA_PBS_GPU_QUEUE:-}" ] && queue="$PMA_PBS_GPU_QUEUE"
    if [ "${PMA_GPU_PROVIDER:-}" = "container" ]; then
        env_opts+=(-v "PMA_DEVICE=gpu")
    else
        env_opts+=(-v "PMA_DEVICE=gpu${PMA_CONDA_ENV_GPU:+,PMA_CONDA_ENV=${PMA_CONDA_ENV_GPU}}")
    fi
else
    # Non-GPU rule: force CPU dispatch when the head job was configured device=gpu.
    env_opts+=(-v "PMA_DEVICE=cpu")
fi

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
[ "${#env_opts[@]}" -gt 0 ] && opts+=("${env_opts[@]}")

# GPU container bind contract — KEEP IDENTICAL in slurm-submit.sh:
#   A container that runs pipeline code MUST bind BOTH the resolved repo root
#   (PMA_REPO_ROOT — for launch_runner.sh + the `executor` package on PYTHONPATH)
#   AND the resolved run directory (PMA_RUN_DIR — for internal/artifacts I/O),
#   because the run data may sit under a nested mount that singularity's default
#   $HOME/$PWD auto-mount does not cover. Optional extra binds (PMA_GPU_BIND,
#   sourced from site.config common.scratch via hpc.env) append after. Producers
#   resolve these paths (launch_runner.sh / configure-execution), so no re-resolve here.
# For a GPU container job, submit a wrapper that exec's the jobscript inside the
# image with GPU passthrough; otherwise submit the jobscript directly.
target="$jobscript"
if [ "${is_gpu}" = "1" ] && [ "${PMA_GPU_PROVIDER:-}" = "container" ] && [ -n "${PMA_GPU_IMAGE:-}" ]; then
    if [ -z "${PMA_RUN_DIR:-}" ]; then
        echo "pbs-submit.sh: WARNING — GPU container job but PMA_RUN_DIR is unset; the run" \
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
    echo "DRY_RUN qsub ${opts[*]} ${target}"
    if [ "$target" != "$jobscript" ]; then
        echo "--- container wrapper ($target) ---"
        cat "$target"
    fi
    exit 0
fi

# qsub -terse outputs just the job id (e.g. "1234567.pbs"); snakemake reads it as
# the submission handle.
exec qsub "${opts[@]}" "$target"
