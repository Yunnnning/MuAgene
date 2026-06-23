# Processing-MuAgent — tool contracts

The `executor` CLI is the **only** state-mutating interface (hard rule 3). Each command
below lists purpose · what it mutates · failure/idempotency. All take `--config <run.yaml>`
(canonical: `deliverables/plan/config/run.yaml` after `init`). State-file lifecycle:
[`../../contracts/state_model.md`](../../contracts/state_model.md). A tripwire test asserts
every command here matches the live CLI.

### executor init
Scaffold the run directory from a draft `run.yaml`; write the canonical config + a blank
`biological_context.md`. Mutates: `internal/`, `deliverables/plan/config/`. Idempotent on an
existing run (won't clobber recorded state).

### executor declare-branch
Record `workflow_branch` (`paired|separate|rna_only|atac_only`) in `parameters.yaml`. Always
confirm the inferred branch with the user first. Mutates: `parameters.yaml`.

### executor configure-execution
Record the execution mode + write `site.config` (+ `hpc.env` for HPC). Requires
`--mode <local|slurm> --confirmed-by-user` (never pass the flag unless the user actually
chose). `--device gpu` only with a cluster mode. Mutates: `parameters.yaml`, `site.config`,
`hpc.env`. Failure: missing `gpu_image_uri` for SLURM GPU → fail loud.

### executor hpc-info
Probe the cluster (partitions/accounts/GPU). **Read-only** — no run state touched.

### executor plan-review
Assemble `preprocessing_plan.json` (defaults from `executor/defaults.py`) and render the
plan-review deliverables + `stage_meta/`. `--intro-context` emits metadata; `--intro "<text>"`
re-renders. Mutates: plan artifacts, `deliverables/plan/plan_review_<run>.md` + `.html`,
`stage_meta/`. Idempotent (deterministic re-render).

### executor approve
Write `internal/checkpoints/<gate>.approved`, unblocking downstream Snakemake rules. Gates:
`plan_review`, `post_qc_review`. Marker-gene flags: `--defer-marker-genes` / `--skip-marker-genes`.
Failure: refuses `plan_review` while the marker-gene decision is unresolved.

### executor revise
Change a planned/QC parameter: `revise <stage> <key>=<value> [--rationale STR]`. Mutates:
`parameters.yaml` (adds `revision_of`) and **deletes** the revised stage's downstream
artifacts so they re-run. At `post_qc_review` this is destructive — diagnose the binding
constraint and confirm first ([`skills/qc_review_and_revise.md`](skills/qc_review_and_revise.md)).
Idempotent: re-revising to the same value is a no-op.

### executor run
Local-only Snakemake. Refuses to launch until `execution.user_confirmed=true`; refuses cluster
modes (that's `submit`). Mutates: stage artifacts + checkpoints as the DAG advances.

### executor submit
Cluster-only. Hands the head-job spec to **Execution-MuAgent** (which submits + supervises);
source `hpc.env` first. Refuses until `execution.user_confirmed=true`. Mutates: `stage_meta/head_job.yaml`,
starts the supervision daemon. After it returns, follow [`skills/hpc_monitoring.md`](skills/hpc_monitoring.md).

### executor status
Per-stage state report (awaiting / running / approved / complete). **Read-only.**

### executor hpc-status
One-shot HPC health, read from Execution's `latest_snapshot.json`. **Read-only.** Never run
`--watch` (blocking) — use report-and-repoll ([`skills/hpc_monitoring.md`](skills/hpc_monitoring.md)).

### executor marker-gene-check
`marker-gene-check <gene1> <gene2> ...`: plot marker expression before/after ambient correction
and refresh the QC reports. `--plot-only` skips the report refresh. Never supply gene names
yourself. Mutates: figures + QC report deliverables.

### executor propose
Regenerate the cheap stage proposal YAMLs (the `*_propose` localrule outputs) without running
heavy compute. Mutates: `internal/proposals/`.

### executor supervisor-restart
Restart the Execution-MuAgent supervision daemon for a live submission **without** resubmitting
the job (recovery from daemon death). Mutates: `internal/hpc_monitor/monitor.pid`.

### executor unlock
Remove a stale Snakemake working-directory lock (after an interrupted run). Mutates:
`internal/snakemake/.snakemake/locks/`.

### executor regenerate-locks
Dev command: re-solve the CPU conda-lock from `workflow/envs/processing.yaml` and stamp its
`# source-sha256:`. Run + commit after editing the env YAML, or `submit`/`validate-env` fail
loud (`lock_stale_vs_yaml`). Mutates: `workflow/envs/processing.linux-64.lock`.
