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

Watch and cancel on detected hang/failure:

```bash
Execution-MuAgent monitor --run-dir /path/to/run --job-id 123456 \
  --watch --stale-minutes 20 --interval 60 --kill-on-hang
```

Reports are written under:

```text
<run_dir>/internal/hpc_monitor/
```

`Processing-MuAgent submit` auto-registers submitted PBS/SLURM head jobs when
this sibling package is present. Set `PMA_HPC_MONITOR_KILL_ON_HANG=1` if the
auto-started watcher should cancel jobs after detecting a stale run.
