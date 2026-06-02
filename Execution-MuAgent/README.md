# Execution-MuAgent

Execution agent for MuAgene. Owns everything between a spec and a running job: renders scheduler submission scripts from the head-job spec and `site.config`, submits to PBS/SLURM, monitors progress using a two-clock state machine, and reports findings back to Processing-MuAgent — without user interaction.

## Architecture

Processing-MuAgent and Execution-MuAgent share a two-file contract:

| File | Written by | Read by | Contains |
|------|-----------|---------|----------|
| `deliverables/pre_run/config/site.config` | Processing-MuAgent | Execution-MuAgent | Platform description: scheduler, partition/queue, account/QOS, conda env or container, resource scale, fs_hang_policy |
| `internal/stage_meta/head_job.yaml` | Processing-MuAgent (`submit`) | Execution-MuAgent | Head-job submission spec: resources (CPU/mem/walltime), input config path, progress_timeout_hint, snakemake_target |
| `internal/stage_meta/<stage>.yaml` | Processing-MuAgent (`plan-review`) | Execution-MuAgent (monitoring) | Per-stage metadata: science_description, resources, inputs/outputs, progress_timeout_hint. Not submitted — used for output validation and monitoring hints. |

`hpc.env` is generated from `site.config` by Processing-MuAgent; it is a shell-variable projection, not an independent source. Do not edit it directly.

Processing-MuAgent never submits jobs directly. `Processing-MuAgent submit` delegates to `Execution-MuAgent execute-spec`, which handles rendering, submission, recording, and monitoring. Snakemake submits per-stage child jobs from within the running head-job.

## Commands

### `execute-spec` — full lifecycle (primary path)

Takes the head-job spec + `site.config`, validates, renders a submission script, submits, records to the execution manifest, registers for monitoring, and optionally watches until the job exits.

```bash
Execution-MuAgent execute-spec \
  --spec /path/to/run/internal/stage_meta/head_job.yaml \
  --site-config /path/to/run/deliverables/pre_run/config/site.config \
  --run-dir /path/to/run \
  --repo-root /path/to/MuAgene/Processing-MuAgent \
  --target all \
  [--watch] [--interval 900] [--kill-on-hang]
```

Steps performed in order:
1. **Validate** — checks resources > 0, scheduler supported, input files exist. On error: writes `spec_validation_error` finding to `latest_report.md` and exits non-zero.
2. **Render** — maps spec resources to scheduler directives (partition, account, QOS, CPU, memory, walltime); wraps command in container invocation if `site.config` specifies one. Writes script to `internal/hpc_monitor/scripts/<stage>_<timestamp>.sh`.
3. **Submit** — `sbatch --parsable` (SLURM) or `qsub -terse` (PBS).
   - **Policy rejection** (invalid partition/account, walltime over site limit): writes `submit_rejected_policy` finding to `latest_report.md`; exits non-zero. Processing-MuAgent relays this as an adjustable hint to the user.
   - **Transient failure**: retries up to 2× with 10 s backoff; reports `submit_rejected_transient` if still failing.
4. **Record** — appends to `internal/hpc_monitor/execution_manifest.jsonl` (stage, science_description, job_id, spec_path, script_path, expected_outputs).
5. **Register** — writes to `internal/hpc_monitor/submissions.jsonl` with `spec_path` and `progress_timeout_hint`.
6. **Monitor** (with `--watch`) — runs the state machine until all jobs exit.

### `register` — break-glass: record a manually-submitted job

Use this when you submitted the head-job manually with `sbatch`/`qsub` (instead of via `Processing-MuAgent submit`). Execution-MuAgent will pick up the job and monitor it as if it had submitted it.

```bash
Execution-MuAgent register \
  --agent Processing-MuAgent \
  --executor slurm \
  --job-id 123456 \
  --run-dir /path/to/run \
  --config /path/to/run/deliverables/pre_run/config/run.yaml \
  --target all \
  --repo-root /path/to/MuAgene/Processing-MuAgent \
  --log-path /path/to/logs/pma_runner-123456.out \
  [--spec-path /path/to/run/internal/stage_meta/head_job.yaml]
```

Pass `--spec-path` so the monitor uses the spec's `progress_timeout_hint` instead of the global `--stale-minutes` fallback.

### `monitor` — watch a registered job

```bash
# Watch with per-spec timeout (preferred):
Execution-MuAgent monitor --run-dir /path/to/run --job-id 123456 \
  --watch --interval 900 \
  --spec-path /path/to/run/internal/stage_meta/head_job.yaml

# Watch with global fallback timeout (--stale-minutes 90 is the default):
Execution-MuAgent monitor --run-dir /path/to/run --job-id 123456 \
  --watch --stale-minutes 90 --interval 900

# One-shot check (no watch):
Execution-MuAgent monitor --run-dir /path/to/run --job-id 123456

# Report only, no cancellation:
Execution-MuAgent monitor --run-dir /path/to/run --job-id 123456 \
  --watch --no-kill-on-hang

# Filesystem-hang kills instead of holding:
Execution-MuAgent monitor --run-dir /path/to/run --job-id 123456 \
  --watch --fs-hang-policy kill_and_resubmit
```

## Monitoring state machine

Detection and decision are always separate. A stall signal is a suspicion, never a verdict.

### Two clocks

**Check interval** (`--interval`, default 900 s / 15 min) — how often the watcher wakes. A sampling rate, the same for every stage. A coarse interval only delays noticing a stall by up to one interval — it never causes a bad kill.

**tolerance_n** — how many consecutive quiet intervals are allowed before raising a stall flag. Derived from the stage's `progress_timeout_hint`: `tolerance_n = ceil(progress_timeout_hint_min × 60 / interval_s)`. The stage declares its tolerance; the interval is just how it is counted.

`progress_timeout_hint` values in `internal/stage_meta/<stage>.yaml` come from `workflow/resources.smk` (the single source of truth), written at `plan-review` time. When no hint is present (e.g. for a manually-registered job without a spec), `--stale-minutes 90` is the fallback default.

A **heartbeat** fires when any run-scoped file mtime advances OR the head log grows since the previous check. Silence resets to 0 on a heartbeat; increments by 1 on a quiet interval.

### States

| State | Meaning | Transition |
|---|---|---|
| `HEALTHY` | No stall signal | → SUSPECT when silence_intervals ≥ tolerance_n |
| `SUSPECT` | Stall flag raised | → INVESTIGATING immediately (same check) |
| `INVESTIGATING` | Gathering evidence | → RECOVERED / CONFIRMED_DEAD / FS_HANG |
| `RECOVERED` | Investigation found life; silence reset | → HEALTHY, continue |
| `CONFIRMED_DEAD` | Evidence confirmed dead | → KILLED (if --kill-on-hang) |
| `FS_HANG` | Filesystem-related hang | → hold or kill per fs_hang_policy |
| `KILLED` | Cancellation sent | → wait for terminal scheduler state |

Definitive signals (`scheduler_failed`, `workflow_error_marker`) bypass the silence counter and go directly to CONFIRMED_DEAD.

### Investigation evidence

Gathered when entering SUSPECT:

| Signal | Source | Alive if… |
|---|---|---|
| Scheduler state | sacct/squeue (already in snapshot) | Not in FAILED_STATES |
| Error markers | Log tails (already in snapshot) | None found |
| CPU utilization | `sstat` (SLURM, best-effort) | CPU > 1 s |
| Memory | `sstat` MaxRSS | RSS > 100 MB |
| Filesystem responsiveness | `os.stat()` in thread, 5 s timeout | Returns within timeout |
| Child storage hang | Child log "Storing output in storage." | Not present |

### Classification rules (first match wins)

1. `scheduler_state` in FAILED_STATES → `confirmed_dead`
2. Error markers in logs → `confirmed_dead`
3. Child RUNNING + "Storing output in storage." → `fs_hang`
4. Filesystem probe timed out → `fs_hang`
5. CPU > 1 s (sstat) → `recovered`
6. Memory > 100 MB (sstat) → `recovered`
7. Filesystem responsive + scheduler RUNNING + no CPU/memory activity → `confirmed_dead`
8. Inconclusive → `recovered` (conservative default)

### Completion validation (D6)

When the head-job reaches terminal COMPLETED state and no problems were detected during the run, Execution-MuAgent validates output artifacts. It reads every `internal/stage_meta/<stage>.yaml` (excluding `head_job.yaml`) and checks that each declared output file exists and is non-empty. Missing or empty outputs are reported as `output_missing` in `latest_report.md` — a COMPLETED scheduler state with missing outputs is treated as a failure.

### Filesystem hang policy

`fs_hang_policy` in site.config (default `"hold"`):
- `"hold"` — writes `filesystem_hang_suspected` finding, does NOT kill. Processing-MuAgent escalates to user. Repeat on subsequent checks until resolved.
- `"kill_and_resubmit"` — kills (children first, head second), routes through normal failure path.

### Finding codes

| Code | Severity | Meaning |
|---|---|---|
| `submit_rejected_policy` | error | Scheduler policy rejection (partition/account/walltime) |
| `scheduler_failed` | error | Scheduler state in terminal failure set |
| `workflow_error_marker` | error | Error keywords in logs |
| `output_missing` | error | Stage output missing or empty after COMPLETED |
| `stall_suspected` | warning | silence_intervals ≥ tolerance_n; entering investigation |
| `stall_confirmed` | error | Investigation confirmed dead — kill was sent |
| `stall_recovered` | warning | Investigation found life; monitoring continues |
| `filesystem_hang_suspected` | error | D-state / degraded storage hang detected |
| `no_progress_files` | warning | No progress files found yet (early stage) |
| `scheduler_completing` | warning | Job in COMPLETING; may indicate epilog/NFS stall |
| `scheduler_query_failed` | warning | Scheduler query timed out or failed |

## Output files

All output lives under `<run_dir>/internal/hpc_monitor/`:

```
internal/hpc_monitor/
├── submissions.jsonl            ← append-only registration log
├── latest_submission.json       ← most recent submission record
├── execution_manifest.jsonl     ← append-only per-submit record (execute-spec path)
├── latest_report.md             ← most recent findings (Processing-MuAgent reads this)
├── latest_snapshot.json         ← full snapshot + monitor_state at last check
├── scripts/
│   └── <stage>_<timestamp>.sh   ← rendered submission scripts
└── reports/
    └── <job_id>_<timestamp>.md  ← historical problem reports
```

`latest_snapshot.json` includes a `"monitor_state"` key with `health`, `silence_intervals`, `tolerance_n`, `investigation`, and `confirmed_dead_reason` — readable by `Processing-MuAgent hpc-status`.

## Head-job spec format

The head-job spec is written by `Processing-MuAgent submit` to `internal/stage_meta/head_job.yaml`:

```yaml
schema_version: '1'
stage: head_job
science_description: Snakemake orchestrator — submits and monitors all per-stage child jobs
resources:
  cpus: 1
  mem_mb: 4000
  walltime_min: 1440
inputs:
  config: /path/to/run/deliverables/pre_run/config/run.yaml
outputs: {}
progress_timeout_hint: 120
snakemake_target: all
```

## Per-stage metadata format

Stage metadata files in `internal/stage_meta/<stage>.yaml` are written by `plan-review`. They are monitoring metadata — not submitted. `progress_timeout_hint` values come from `workflow/resources.smk`:

```yaml
schema_version: '1'
stage: s3_doublets
science_description: Detect and remove doublets using Scrublet (RNA) and SnapATAC2 (ATAC)
resources:
  cpus: 2
  mem_mb: 32000
  walltime_min: 330
inputs:
  rna_h5ad: /path/to/run/internal/artifacts/s1_rna_qc/rna_qc.h5ad
  atac_h5ad: /path/to/run/internal/artifacts/s2_atac_qc/atac_qc.h5ad
outputs:
  rna_post: /path/to/run/internal/artifacts/s3_doublets/rna_post_doublet.h5ad
  atac_post: /path/to/run/internal/artifacts/s3_doublets/atac_post_doublet.h5ad
progress_timeout_hint: 60
```

## site.config format

`site.config` is the single platform source of truth. Written by `Processing-MuAgent configure-execution`:

```yaml
schema_version: '1'
scheduler: slurm       # slurm | pbs
slurm:
  partition: cpu-medium
  account: vaquerizas-lab
  qos: null
common:
  resources_scale: 2
  conda_env: muagene
  container: null      # path to .sif for Apptainer/Singularity, or null
  scratch: null
  fs_hang_policy: hold # hold | kill_and_resubmit
```

`hpc.env` is derived from this file by Processing-MuAgent and cannot drift from it.

## Install

```bash
cd /path/to/Execution-MuAgent
pip install -e .
```

Requires Python 3.10+, `click>=8.1`, `pyyaml>=6`.
