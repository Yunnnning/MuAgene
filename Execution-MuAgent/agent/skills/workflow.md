# Execution-MuAgent — operational workflow

This is the canonical procedure for runtime execution and machine provisioning. Finding
codes and state ownership are defined in [`../../../contracts/`](../../../contracts/);
command mutations and flags are defined in [`../tools.md`](../tools.md).

## `execute-spec` — Processing-invoked runtime

Processing calls `execute-spec` after writing `site.config` and the head/stage specs.
Execution never changes those inputs.

### 1. Validate and preflight

1. Load the head-job spec and `site.config`.
2. Validate positive resources, SLURM support, and declared inputs.
3. Reconcile the CPU environment. Reconcile GPU only when the site selects GPU and the
   spec declares GPU-capable stages.
4. If validation or environment preflight fails:
   - preserve the original registered environment finding when available;
   - otherwise use the appropriate validation/provisioning finding;
   - write it to `latest_snapshot.json`;
   - exit non-zero without rendering or submitting.

Never convert an environment error into a warning or fall back to a different device.

### 2. Render and submit

1. Map the spec resources and confirmed site settings to SLURM directives.
2. Wrap the command with the configured container invocation when required.
3. Write the scheduler script and call `sbatch --parsable`.
4. Classify rejection before retrying:
   - policy rejection: write `submit_rejected_policy` to the snapshot and stop;
   - transient rejection: retry at most twice with the defined backoff, then write
     `submit_rejected_transient` and stop.
5. On success, record the job ID, script, expected outputs, and submission metadata; update
   the latest-submission record used by `resume-monitor`.

The optional Markdown report is debug/audit output. Processing consumes only the structured
snapshot.

### 3. Monitor when `--watch`

Use the stage's `progress_timeout_hint` and the monitor's canonical check interval.
Refresh `latest_snapshot.json` on every check, including healthy checks.

For each check:

1. Collect bounded SLURM (`sacct`/`squeue`) state, run-scoped file/log progress, child-job
   state, and error markers.
2. Verify any declared outputs that have appeared and emit one informational verification
   finding per completed stage.
3. Treat scheduler failure and workflow error markers as definitive unhealthy evidence.
4. If the workflow is cleanly complete but the head process lingers, cancel the head only,
   emit completion, and exit without a failure action.
5. Otherwise, count quiet intervals. Once the stage-specific tolerance is reached,
   investigate; silence alone never authorizes cancellation.
6. During investigation, classify the first conclusive signal:
   - failed scheduler state or error marker → `CONFIRMED_DEAD`;
   - storage sentinel or unresponsive filesystem → `FS_HANG`;
   - advancing CPU time, or GPU utilization when available → `RECOVERED`;
   - responsive filesystem plus RUNNING plus a measured-flat CPU sample →
     `CONFIRMED_DEAD`;
   - missing first comparison sample → `RECOVERED` and establish a baseline.
7. For `CONFIRMED_DEAD` or `FS_HANG`, cancel children before the head, record the reason,
   persist the snapshot, and exit. Never hold for a human or resubmit.

Processing creates `monitor.pid` before launching the daemon. Remove it in `finally` on
normal completion, failure, or cancellation.

## `resume-monitor` — Processing-invoked recovery

Use only when the supervisor stopped but the recorded cluster job is still active.

1. Load the latest successful submission record.
2. Resume the same watch loop without submitting another job.
3. Remove the Processing-created PID file in `finally`.

## Machine provisioning — operator-facing

These commands may print directly to the operator because they run before a Processing
session exists.

### `init-machine`

1. Probe the Linux host for a conda-compatible manager, SLURM, container runtime, and GPU.
2. Write the machine profile.
3. Create/update the integrated `muagene` CPU environment from the committed lock.
4. Install both agent packages into that environment without re-resolving dependencies.
5. When requested, pull the pinned GPU image; never build it locally.
6. Validate requested devices and print a readiness report. Fail loud on any error.

### `provision-env` / `validate-env` / `doctor`

- Provisioning is idempotent and recipe-driven: CPU uses the committed lock; GPU pulls a
  pinned image.
- Validation checks presence, fingerprints, and declared imports without mutating the env.
- `doctor` reports capabilities and validation results.
- A run-specific site configuration is optional after machine bootstrap; otherwise use the
  machine profile plus the shared environment manifest.

Do not duplicate failure-code explanations here. Emit the registered code and actionable
message from the canonical finding registry.
