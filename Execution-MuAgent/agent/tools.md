# Execution-MuAgent — tool contracts

Each command lists purpose · what it mutates · failure/idempotency. Run-time state is written
under `internal/hpc_monitor/`; environment reconciliation may also update
`~/.muagene/env_state.json`. Bootstrap commands are operator-facing and print to stdout.
Finding codes: [`../../contracts/findings.yaml`](../../contracts/findings.yaml);
state lifecycle: [`../../contracts/state_model.md`](../../contracts/state_model.md). A tripwire
test asserts every command here matches the live CLI.

### Execution-MuAgent execute-spec
The run-time lifecycle entry point (Processing invokes it): validate spec → render submission
script → submit (`sbatch`) → record, then optionally monitor with `--watch`. Processing supplies
the spec, site config, run directory, repository root, and target; monitoring controls are
`--watch`, `--interval`, and `--kill-on-hang`. Every classified pre-submit failure writes a
snapshot finding and exits non-zero; transient submission retries ≤2×. Never resubmits.

### Execution-MuAgent resume-monitor
Resume supervision of the latest submission **without** resubmitting (recovers from daemon
death). Reads `latest_submission.json` and re-enters the watch loop. Processing creates the
supervisor PID file; this command removes it on exit.

### Execution-MuAgent report
Print the latest `latest_report.md` for manual debugging. It is read-only and is not the
Processing contract; Processing reads `latest_snapshot.json`.

### Execution-MuAgent init-machine
Fresh-machine provisioning (operator-facing): probe capabilities → write
`~/.muagene/machine.config` → reconcile the requested environment → install both packages →
pull a GPU image when GPU setup is requested → validate → print readiness. Mutates:
`~/.muagene/` and the configured environment.

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
Print machine capabilities (scheduler, env manager, container runtime, GPU presence). If a
machine profile exists, also validate its configured environments. **Read-only** diagnostics.
