# Execution-MuAgent (Internal Runtime Orchestrator)

The Execution Agent operates entirely in the background as MuAgene's execution layer. It **translates** workflow specifications into scheduler-ready jobs, **submits** them to PBS/SLURM, and **supervises** execution using a dual-clock state machine (heartbeat signals + walltime tracking), and **reports** the job status to Processing-MuAgent. By continuously collecting telemetry, identifying failures, detecting stalled jobs, and generating diagnostics, it maintains execution observability and traceability throughout the job lifecycle. It is an internal support agent with **no direct user interaction**.

```
Head-job Spec → Job Rendering → Cluster Submission → Runtime Monitoring → Diagnostics → Status Report.
```

## Architecture

Processing-MuAgent and Execution-MuAgent share a two-file contract:

| File | Written by | Read by | Contains |
|------|-----------|---------|----------|
| `deliverables/plan/config/site.config` | Processing-MuAgent | Execution-MuAgent | Platform description: scheduler, partition/queue, account/QOS, conda env or container, resource scale |
| `internal/stage_meta/head_job.yaml` | Processing-MuAgent (`submit`) | Execution-MuAgent | Head-job submission spec: resources (CPU/mem/walltime), input config path, progress_timeout_hint, snakemake_target |
| `internal/stage_meta/<stage>.yaml` | Processing-MuAgent (`plan-review`) | Execution-MuAgent (monitoring) | Per-stage metadata: science_description, resources, inputs/outputs, progress_timeout_hint. Not submitted — used for output validation and monitoring hints. |

`hpc.env` is generated from `site.config` by Processing-MuAgent; it is a shell-variable projection, not an independent source. Do not edit it directly.

Processing-MuAgent never submits or monitors cluster jobs directly — `Processing-MuAgent run` is local-only and `Processing-MuAgent submit` is cluster-only. `submit` delegates to `Execution-MuAgent execute-spec`, which handles rendering, submission, recording, and monitoring. Snakemake submits per-stage child jobs from within the running head-job.

## Commands

### `execute-spec` — full lifecycle (primary path)

Takes the head-job spec + `site.config`, validates, renders a submission script, submits, records to the execution manifest, registers for monitoring, and optionally watches until the job exits.

```bash
Execution-MuAgent execute-spec \
  --spec /path/to/run/internal/stage_meta/head_job.yaml \
  --site-config /path/to/run/deliverables/plan/config/site.config \
  --run-dir /path/to/run \
  --repo-root /path/to/MuAgene/Processing-MuAgent \
  --target all \
  [--watch] [--interval 270] [--kill-on-hang]
```

Steps performed in order:
1. **Validate** — checks resources > 0, scheduler supported, input files exist. On error: writes `spec_validation_error` finding to `latest_report.md` and exits non-zero.
2. **Render** — maps spec resources to scheduler directives (partition, account, QOS, CPU, memory, walltime); wraps command in container invocation if `site.config` specifies one. Writes script to `internal/hpc_monitor/scripts/<stage>_<timestamp>.sh`.
3. **Submit** — `sbatch --parsable` (SLURM) or `qsub -terse` (PBS).
   - **Policy rejection** (invalid partition/account, walltime over site limit): writes `submit_rejected_policy` finding to `latest_report.md`; exits non-zero. Processing-MuAgent relays this as an adjustable hint to the user.
   - **Transient failure**: retries up to 2× with 10 s backoff; reports `submit_rejected_transient` if still failing.
4. **Record** — appends to `internal/hpc_monitor/execution_manifest.jsonl` (stage, science_description, job_id, spec_path, script_path, expected_outputs).
5. **Register** — writes to `internal/hpc_monitor/submissions.jsonl` with `spec_path` and `progress_timeout_hint`.
6. **Monitor** (with `--watch`) — runs the watch loop until all jobs exit, then removes `internal/hpc_monitor/monitor.pid`.

### `resume-monitor` — restart supervision without resubmitting

Reads `internal/hpc_monitor/latest_submission.json`, reconstructs the monitoring context, and runs the same watch loop as `execute-spec --watch` — but without submitting a new job. Invoked by `Processing-MuAgent supervisor-restart` when the supervision daemon dies mid-run.

```bash
Execution-MuAgent resume-monitor \
  --run-dir /path/to/run \
  [--interval 270] [--kill-on-hang]
```

This is not a command you call directly — Processing-MuAgent starts it as a background daemon. It removes `monitor.pid` when it finishes (whether the job completed, was killed, or the command crashed).

### `report` — read the latest diagnostic report

Prints `<run_dir>/internal/hpc_monitor/latest_report.md` — the findings and diagnostics from the most recent monitor check.

```bash
Execution-MuAgent report --run-dir /path/to/run
```

`report` is the only command a human uses to inspect what Execution-MuAgent found. `execute-spec` and `resume-monitor` are machine entry points invoked by Processing-MuAgent. There are no manual-submission, registration, or standalone-monitor commands beyond these three.

## Monitoring state machine

Detection and decision are always separate. A stall signal is a suspicion, never a verdict.

### Two clocks

**Check interval** (`--interval`, default 270 s / 4.5 min — the constant `monitor.DEFAULT_CHECK_INTERVAL_S`) — how often the watcher wakes. A sampling rate, the same for every stage. A coarse interval only delays noticing a stall by up to one interval — it never causes a bad kill. Every snapshot records `interval_s` and `next_recheck_after_s` (= interval + `REPOLL_BUFFER_S`, ~25 s) so Processing-MuAgent re-polls just after each check rather than hardcoding a cadence.

**tolerance_n** — how many consecutive quiet intervals are allowed before raising a stall flag. Derived from the stage's `progress_timeout_hint`: `tolerance_n = ceil(progress_timeout_hint_min × 60 / interval_s)`. The stage declares its tolerance; the interval is just how it is counted.

`progress_timeout_hint` values in `internal/stage_meta/<stage>.yaml` come from `workflow/resources.smk` (the single source of truth), written at `plan-review` time. When no hint is present (e.g. for the head-job spec itself), a 90-minute fallback is used.

A **heartbeat** fires when any run-scoped file mtime advances OR the head log grows since the previous check. Silence resets to 0 on a heartbeat; increments by 1 on a quiet interval.

### States

| State | Meaning | Transition |
|---|---|---|
| `HEALTHY` | No stall signal | → SUSPECT when silence_intervals ≥ tolerance_n |
| `SUSPECT` | Stall flag raised | → INVESTIGATING immediately (same check) |
| `INVESTIGATING` | Gathering evidence | → RECOVERED / CONFIRMED_DEAD / FS_HANG |
| `RECOVERED` | Investigation found life; silence reset | → HEALTHY, continue |
| `CONFIRMED_DEAD` | Evidence confirmed dead | → KILLED (if --kill-on-hang) |
| `FS_HANG` | Filesystem-related hang | → KILLED (if --kill-on-hang); reported to Processing |
| `KILLED` | Cancellation sent | → wait for terminal scheduler state |

Definitive signals (`scheduler_failed`, `workflow_error_marker`) bypass the silence counter and go directly to CONFIRMED_DEAD.

An unhealthy verdict (`CONFIRMED_DEAD` or `FS_HANG`) is killed for cleanup and reported. Execution-MuAgent never holds for a human and never resubmits — Processing-MuAgent reads the report, escalates to the human, fixes, and resubmits.

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

### Output verification (per step + terminal)

Output verification is proper — not a folder/size check. `verify_output_file` opens each declared output and confirms it is a complete, loadable file:
- `.h5ad`/`.h5mu`/`.h5` — opened with `h5py` and checked for the expected root groups (`X`/`obs`/`var`, or `mod` for h5mu). Without `h5py`, the HDF5 8-byte signature is checked at superblock offsets.
- `.parquet` — `pyarrow.parquet.read_metadata` (footer); fallback to the `PAR1` head/tail magic.
- `.json` — parsed; text sentinels — non-empty.

Because stages write outputs atomically (`/tmp` stage + fsync + `os.rename`), a file present at its final path is complete, so a valid signature plus non-zero size is a strong correctness signal.

**Per step (normal progress):** on every check the monitor verifies each spec's declared outputs that have appeared and emits a one-time `stage_output_verified` finding. The head-job spec is verified like any other when Processing has populated its `outputs` from the target stage (e.g. a planning `s0_ingest_execute` submission), so a clean head-job exit still emits this finding instead of leaving `verified_stages` empty. `latest_snapshot.json` (with `monitor_state.verified_stages`) is refreshed every check — healthy or not — so Processing's `hpc-status` never reads stale state.

**Terminal:** when the head-job reaches COMPLETED, `validate_terminal_outputs` runs the same verifier over every `internal/stage_meta/<stage>.yaml` (excluding `head_job.yaml`). Any missing, empty, or corrupt output is reported as `output_missing` — a COMPLETED scheduler state with an unverifiable output is treated as a failure.

### Unhealthy runs: no human fallback in Execution

When a run is unhealthy (`CONFIRMED_DEAD` or `FS_HANG`), Execution-MuAgent kills the job (children first, then head) for cleanup and writes diagnostics. It never holds for a human and never resubmits. Processing-MuAgent reads the report, reports to the human, implements the fix, and resubmits. There is no `fs_hang_policy` knob — filesystem hangs follow the same kill-and-report path as any other confirmed-dead verdict.

### Finding codes

| Code | Severity | Meaning |
|---|---|---|
| `submit_rejected_policy` | error | Scheduler policy rejection (partition/account/walltime) |
| `scheduler_failed` | error | Scheduler state in terminal failure set |
| `workflow_error_marker` | error | Error keywords in logs |
| `output_missing` | error | Stage output missing, empty, or corrupt after COMPLETED |
| `stage_output_verified` | info | A stage's declared outputs verified complete and loadable |
| `stall_suspected` | warning | silence_intervals ≥ tolerance_n; entering investigation |
| `stall_confirmed` | error | Investigation confirmed dead — kill was sent |
| `stall_recovered` | warning | Investigation found life; monitoring continues |
| `filesystem_hang_suspected` | error | D-state / degraded storage hang — kill was sent |
| `no_progress_files` | warning | No progress files found yet (early stage) |
| `scheduler_completing` | warning | Job in COMPLETING; may indicate epilog/NFS stall |
| `scheduler_query_failed` | warning | Scheduler query timed out or failed |

## Output files

All output lives under `<run_dir>/internal/hpc_monitor/`:

```
internal/hpc_monitor/
├── submissions.jsonl            ← append-only registration log
├── latest_submission.json       ← most recent submission record (used by resume-monitor)
├── execution_manifest.jsonl     ← append-only per-submit record (execute-spec path)
├── latest_report.md             ← most recent findings (Processing-MuAgent reads this)
├── latest_snapshot.json         ← full snapshot + monitor_state at last check
├── monitor.pid                  ← PID of the running supervision daemon (removed on exit)
├── monitor.log                  ← symlink → most recent monitor_<timestamp>.log
├── monitor_<timestamp>.log      ← daemon output for each submit or supervisor-restart
├── scripts/
│   └── <stage>_<timestamp>.sh   ← rendered submission scripts
└── reports/
    └── <job_id>_<timestamp>.md  ← historical problem reports
```

`latest_snapshot.json` includes a `"monitor_state"` key with `health`, `silence_intervals`, `tolerance_n`, `investigation`, `confirmed_dead_reason`, and `verified_stages` — readable by `Processing-MuAgent hpc-status`. It also carries `interval_s` and `next_recheck_after_s` (the daemon's poll cadence + re-poll buffer) so Processing's re-poll delay is data-driven. It is refreshed on every check, including healthy ones.

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
  config: /path/to/run/deliverables/plan/config/run.yaml
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

S1/S2 declare their `qc_summary.json` as the monitored output (not the large h5ad), because the
h5ad is deleted by `executor approve post_qc_review` and would cause false `output_missing`
findings in subsequent post-QC job runs:

```yaml
# s1_rna_qc.yaml
outputs:
  qc_summary_json: /path/to/run/internal/artifacts/s1_rna_qc/qc_summary.json

# s2_atac_qc.yaml
outputs:
  qc_summary_json: /path/to/run/internal/artifacts/s2_atac_qc/qc_summary.json
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
```

`hpc.env` is derived from this file by Processing-MuAgent and cannot drift from it.

## Install

```bash
cd /path/to/Execution-MuAgent
pip install -e .
```

Requires Python 3.10+, `click>=8.1`, `pyyaml>=6`.
