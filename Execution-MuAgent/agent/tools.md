# Execution-MuAgent — tool contracts

Each command lists purpose · what it mutates · failure/idempotency. Run-time commands write
only to `internal/hpc_monitor/` (Processing is the sole consumer); bootstrap commands are
operator-facing and print to stdout. Finding codes: [`../../contracts/findings.yaml`](../../contracts/findings.yaml);
state lifecycle: [`../../contracts/state_model.md`](../../contracts/state_model.md). A tripwire
test asserts every command here matches the live CLI.

### Execution-MuAgent execute-spec
The run-time lifecycle entry point (Processing invokes it): validate spec → render submission
script → submit (`sbatch`/`qsub`) → record → monitor until exit. **No user-facing output** —
all state goes to `internal/hpc_monitor/` (`latest_snapshot.json`, `execution_manifest.jsonl`,
`submissions.jsonl`, `scripts/`, `monitor.pid`). Flags: `--watch`, `--interval`, `--kill-on-hang`.
Failure: policy rejection → finding + non-zero exit; transient → retry ≤2×. Never resubmits.

### Execution-MuAgent resume-monitor
Resume supervision of the latest submission **without** resubmitting (recovers from daemon
death). Reads `latest_submission.json`; re-arms the watch loop; refreshes `monitor.pid`.

### Execution-MuAgent report
Print the latest `latest_report.md` — an **internal debug/audit helper only**. Processing never
parses it and it is never shown to the user. Read-only.

### Execution-MuAgent init-machine
Fresh-machine bootstrap (operator-facing): probe capabilities → write `~/.muagene/machine.config`
→ create the CPU env from the committed conda-lock → `pip install -e` both packages → pull the
GPU image → validate → print a readiness report. Mutates: `~/.muagene/`, the conda env.

### Execution-MuAgent provision-env
`--device cpu|gpu|both`: idempotently create/update the env from the `environments:` recipe (CPU
conda-lock; GPU pull-only container) and record a fingerprint in `~/.muagene/env_state.json`.
`--site-config` optional once `machine.config` exists. Failure: `platform_unsupported` (non-linux
CPU lock), `gpu_image_unavailable`, `provision_failed` — all fail loud, never degrade.

### Execution-MuAgent validate-env
`--device cpu|gpu|both`: confirm the env is present and import-check its declared modules
(GPU: `rapids_singlecell`/`cupy`). **Read-only**; fails loud (`import_failed`, `lock_stale_vs_yaml`,
`env_missing`/`env_stale`). On a CPU login node the GPU import check defers (`gpu_import_needs_node`).

### Execution-MuAgent doctor
Print machine capabilities (scheduler, env manager, container runtime, GPU presence) and validate
the configured envs. **Read-only** diagnostics.
