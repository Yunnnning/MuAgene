# Execution-MuAgent

Central execution monitor for MuAgene subagents.

Subagents register HPC submissions here so one watchdog can detect stalled
cluster runs, write diagnostics, and optionally cancel jobs.

## Commands

Register a submission:

```bash
Execution-MuAgent register \
  --agent Processing-MuAgent \
  --executor slurm \
  --job-id 123456 \
  --run-dir /path/to/run \
  --config /path/to/run/deliverables/pre_run/config/run.yaml \
  --target post_qc_review_propose \
  --repo-root /path/to/MuAgene/Processing-MuAgent \
  --log-path /path/to/MuAgene/Processing-MuAgent/logs/pma_runner-123456.out
```

Monitor once:

```bash
Execution-MuAgent monitor --run-dir /path/to/run --job-id 123456
```

Watch and cancel on detected hang/failure (default):

```bash
Execution-MuAgent monitor --run-dir /path/to/run --job-id 123456 \
  --watch --stale-minutes 20 --interval 60
```

Report only (no scancel):

```bash
Execution-MuAgent monitor --run-dir /path/to/run --job-id 123456 \
  --watch --no-kill-on-hang
```

Reports are written under:

```text
<run_dir>/internal/hpc_monitor/
```

`Processing-MuAgent submit` auto-registers submitted PBS/SLURM head jobs when
this sibling package is present. The auto-started watcher cancels hung jobs by
default (`PMA_HPC_MONITOR_KILL_ON_HANG=0` disables cancellation). Set
`PMA_NOTIFY_EMAIL` to receive hang alerts in addition to runner exit mail.
