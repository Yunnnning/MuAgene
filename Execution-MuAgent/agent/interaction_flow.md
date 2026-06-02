# Execution-MuAgent — interaction flow

Execution-MuAgent is invoked programmatically by Processing-MuAgent. It never interacts with the user. This document describes the full spec → running job lifecycle.

---

## Triggered by: `execute-spec`

Processing-MuAgent calls this after it has written `site.config` and the per-stage spec.

### Step 1 — Load and validate

1. Load `site.config` → `SiteConfig` (scheduler, partition, account, QOS, container, env).
2. Load `internal/specs/<stage>.yaml` → `StageSpec` (resources, inputs, outputs, progress_timeout_hint, science_description).
3. `validate_spec(spec, site_config)`:
   - Resources > 0.
   - Scheduler is supported.
   - Declared input paths exist on the filesystem.
4. On any validation error: write `spec_validation_error` finding to `internal/hpc_monitor/latest_report.md`; exit non-zero. Processing-MuAgent surfaces the error to the user.

### Step 2 — Render submission script

1. Map spec resources → scheduler header directives:
   - SLURM: `#SBATCH --cpus-per-task`, `--mem`, `--time`, `--partition`, `--account`, `--qos`
   - PBS: `#PBS -l select=1:ncpus=...:mem=...mb`, `-l walltime=`, `-q`, `-P`
2. If `site_config.container` is set: wrap command in `apptainer exec --bind <run_dir> <container> bash <launch_runner.sh>`.
3. Write script to `internal/hpc_monitor/scripts/<stage>_<timestamp>.sh`.

### Step 3 — Submit and diagnose

1. Run `sbatch --parsable <script>` or `qsub -terse <script>`.
2. On success: capture job_id from stdout.
3. On failure:
   - **Policy rejection** (invalid partition/account, walltime over limit): write `submit_rejected_policy` finding to `latest_report.md`; exit non-zero. Processing-MuAgent reads this and tells the user what to adjust.
   - **Transient failure**: retry up to 2× with 10 s backoff. If still failing: write `submit_rejected_transient` finding; exit non-zero.

### Step 4 — Record to execution manifest

Append to `internal/hpc_monitor/execution_manifest.jsonl`:

```json
{
  "submitted_at": "2026-06-02T14:00:00Z",
  "stage": "s3_doublets",
  "science_description": "Detect and remove doublets using Scrublet (RNA) and SnapATAC2 (ATAC)",
  "job_id": "987654",
  "spec_path": ".../internal/specs/s3_doublets.yaml",
  "script_path": ".../internal/hpc_monitor/scripts/s3_doublets_20260602T140000Z.sh",
  "run_dir": "/path/to/run",
  "expected_outputs": { "rna_post": "...", "atac_post": "..." }
}
```

### Step 5 — Register for monitoring

Write to `internal/hpc_monitor/submissions.jsonl` with `spec_path` and `progress_timeout_hint` fields so the watch loop uses the per-spec timeout.

### Step 6 — Monitor (if `--watch`)

Poll loop at `interval_s` (default 900 s / 15 min). The monitor drives a state machine; detection and decision are always separate.

**Two clocks:**
- `interval_s` — sampling rate (how often the watcher wakes). Same for every stage.
- `tolerance_n = ceil(progress_timeout_hint_min * 60 / interval_s)` — how many consecutive quiet intervals are allowed. A heartbeat fires when any run-scoped file mtime advances OR the head log size grows.

**State machine per iteration:**

1. **Collect** — `collect_snapshot()`: scheduler state (sacct/squeue/qstat, ≤5 s), filesystem progress files, log tails, child job IDs, error markers.

2. **Definitive signals** (always checked, bypass silence machine):
   - `scheduler_failed` → CONFIRMED_DEAD immediately.
   - `workflow_error_marker` (Traceback, OOM, WorkflowError, …) → CONFIRMED_DEAD immediately.

3. **Watcher** (HEALTHY / RECOVERED state): if latest_progress_file mtime or head_log size grew since last check → heartbeat, silence resets to 0. Otherwise `silence_intervals += 1`. When `silence_intervals >= tolerance_n` → SUSPECT, emit `stall_suspected` (warning only).

4. **Investigation** (SUSPECT state): gather independent evidence — `sstat` CPU/memory (SLURM), filesystem responsiveness probe (D-state detection), child storage-hang sentinel ("Storing output in storage."). Classify by rules:
   - Scheduler failed / error markers → `confirmed_dead`
   - Child storage hang / filesystem probe timeout → `fs_hang`
   - CPU active / memory active → `recovered`
   - All silent + responsive filesystem + RUNNING → `confirmed_dead`
   - Inconclusive → `recovered` (conservative default)

5. **Kill** (CONFIRMED_DEAD only):
   - `cancel_submission_jobs()` — children first, then head.
   - Record `confirmed_dead_reason` in cancel result.
   - Write report to `latest_report.md` and `reports/<job_id>_<timestamp>.md`.
   - Write full snapshot + monitor state to `latest_snapshot.json`.

6. **FS hang** (FS_HANG state):
   - Policy `"hold"` (default): write `filesystem_hang_suspected` finding; do NOT kill; report to Processing-MuAgent; continue monitoring.
   - Policy `"kill_and_resubmit"`: kill and route through normal failure path.

7. **Loop exit** when no jobs are active in scheduler (all terminal states).

---

## Triggered by: `register` + `monitor` (legacy/fallback path)

When Processing-MuAgent did direct `sbatch`/`qsub` submission (Execution-MuAgent unavailable at submit time, or fallback path), it calls `register` + starts `monitor --watch` as a background daemon.

- `register`: records submission metadata to `submissions.jsonl`; accepts optional `--spec-path` to link the spec for per-spec timeout.
- `monitor --watch`: same watch loop as Step 6 above; loads spec timeout via `--spec-path` or from the registered submission record.

---

## Reporting back to Processing-MuAgent

All findings land in `internal/hpc_monitor/latest_report.md`; monitor state in `latest_snapshot.json`. Processing-MuAgent reads these:

| Finding code | Meaning | Processing-MuAgent action |
|---|---|---|
| `submit_rejected_policy` | Scheduler rejected submission (partition/account/walltime) | Relay as adjustable hint; ask user to fix site.config |
| `scheduler_failed` | Scheduler reports terminal failure state | Relay finding; route to failure path |
| `workflow_error_marker` | Error keywords in logs (Traceback, OOM, …) | Relay finding; route to failure path |
| `stall_suspected` | N quiet intervals; investigation starting | Informational; no user action needed |
| `stall_confirmed` | Investigation concluded confirmed dead | Relay reason; route to failure path |
| `stall_recovered` | Investigation found life; monitoring continues | Informational; no action needed |
| `filesystem_hang_suspected` | D-state / storage-degraded hang (policy=hold) | Escalate to user; hold until resolved |
| Clean exit | No findings; job completed normally | Continue normal pipeline flow |
