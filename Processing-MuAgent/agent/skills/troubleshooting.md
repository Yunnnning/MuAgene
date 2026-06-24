---
name: troubleshooting
domain: recovery
purpose: Map an observed error / finding to its stage and the executor remedy. Loaded only when something goes wrong.
activation: a finding/error in the snapshot, or an executor command raised
inputs: [internal/hpc_monitor/latest_snapshot.json, internal/artifacts/s0_ingest/validation_report.json]
outputs: []
calls_tools: [hpc-status, configure-execution, revise, declare-branch, submit, run, plan-review, supervisor-restart]
reads_contracts: [findings, latest_snapshot]
writes_state: [parameters.yaml]
handoff: { next: re-enter index.md router, when: remedy applied, on_error: troubleshooting }
---

# Troubleshooting — symptom → remedy

Loaded on demand when a stage skill hits an error. Apply the remedy, then re-run
`executor status` / `hpc-status` and **re-enter [`index.md`](index.md)**. Never retry a logic
error silently (hard rule 5) — relay the raised message, root-cause, then act. Finding-code
meanings + recovery actions are in [`../../../contracts/findings.yaml`](../../../contracts/findings.yaml).

## Ingest / branch (S0)

- **S0 auto-downgrades `paired → separate`** (not an error). No ladder rung (direct overlap,
  suffix-normalized, `barcode_translation_path`, `cell_metadata_path` with `atac_barcode`)
  validated cell-level pairing. Surface `internal/artifacts/s0_ingest/validation_report.json`
  verbatim (`pairing.downgrade_reason` has the specifics). Ask: (a) proceed on `separate`,
  (b) supply `barcode_translation_path` + rerun S0, or (c) abort and fix inputs upstream.
- **S0 raises "declared=… conflicts with detected=…"** — `declare-branch` disagrees with S0
  detection and the declaration is single-modality (`rna_only`/`atac_only`). Relay the
  message; ask the user to correct the declaration or drop the unwanted modality from config.
- **S0 raises "pairing is ambiguous"** — RNA+ATAC Jaccard 30–80% with no subset relation and
  no resolving declaration. Ask paired vs separate; `executor declare-branch <paired|separate>`
  then re-run; for `paired`, supply `barcode_translation_path` if whitelists differ. Don't auto-pick.
- **S3 raises "paired-branch joint barcode intersection is empty"** — S0 committed `paired`
  but no cell survived both modalities' QC + doublet removal (usually QC too aggressive).
  Ask the user to revise S1/S2 thresholds (`executor revise s1_rna_qc …` / `s2_atac_qc …`).
  If pairing used `pairing.translation_table`, check it covers the QC-surviving cells.
- **Re-processing a previously-approved run (rna_ingest.h5ad was cleaned)** — not an error.
  `rna_ingest.h5ad` is a deletable S0 *cache* the post-QC cleanup removes; S1a reconstructs it
  deterministically from the original input via `io.load_rna_ingest` (logged
  `rna_ingest_reconstructed`), so re-processing needs **no S0 re-run**. Just reset the downstream
  artifacts you want regenerated, `revise` the threshold (which auto-refreshes the previews +
  plan files), and resubmit. (Only the truly missing case — `rna_ingest.h5ad` absent **and**
  run.yaml has no `rna_path` — raises, asking you to re-run S0.)

## Context / execution-mode gates

- **Phase 1 gate "biological_context.md is empty"** — user gave no context and didn't opt
  out. Ask for context OR offer the explicit opt-out `--no-context` on the entry point
  (`executor run … --target plan_review_propose --no-context` local, or
  `executor submit … --no-context` on HPC).
- **`run`/`submit` "Execution mode is not set" / "not confirmed by the user"** — compute was
  launched before confirming local vs HPC (fires on fresh and resumed runs). Confirm the
  mode with the user, probe `executor hpc-info` for clusters, then record their explicit
  choice: `executor configure-execution --config $CFG --mode <local|slurm> --confirmed-by-user`.
  Never pass `--confirmed-by-user` without having asked. Re-run the same command.
- **`run` "execution.mode is 'slurm' but `run` is local-only"** — source `hpc.env` and use
  `executor submit --config $CFG --executor slurm`.

## Resources / scheduler / environment (HPC)

- **S0 OOMs / Killed / hits walltime (HPC)** — resource-sizing, not location. Raise
  `PMA_RESOURCES_SCALE` (`configure-execution --mode slurm --resources-scale N …`), then
  `executor submit --config $CFG --executor slurm` again (omit `--target`). (Re-config of the
  *same* mode preserves the existing confirmation — no `--confirmed-by-user` for a
  resource-only change.) **Local-mode** S0 OOM → the machine is too small; switch to HPC. No
  automatic local→cluster retry.
- **A stage execute fails at runtime** — relay the failure (HPC: one-shot `executor
  hpc-status --config $CFG` renders the daemon's structured findings; local: snakemake
  stderr). Root-cause first; don't retry silently. On user insistence, re-`executor submit
  --executor slurm --target <stage>_execute` (HPC) or `executor run --target <stage>_execute`
  (local).
- **S2 fails with `H5Fcreate(): unable to lock file, errno = 11`** — stale `atac_snap.h5ad`
  from a previously-killed run holds an HDF5 POSIX lock; SnapATAC2 cannot recreate it.
  Fix: `rm internal/artifacts/s2_atac_qc/atac_snap.h5ad` (untracked; safe to delete),
  then `executor submit`. The code now auto-deletes this file at stage start (so this
  only arises on an older code version). Also clear any `snapatac2_tmp/` temp subdirs.
- **`submit_rejected_policy`** — scheduler rejected the job (invalid partition/account, or
  walltime over site limit). One-shot `hpc-status` shows the scheduler's exact message. Fix
  the field: partition/account via `executor configure-execution --mode <scheduler> …`
  (rewrites `site.config`), or walltime by reducing `PMA_RESOURCES_SCALE`. Then `submit` again.
- **Environment-preflight error at submit** — provisioning is owned by Execution-MuAgent;
  `submit` auto-provisions a missing/stale env (policy=auto) but fails loud rather than
  degrade. Relay the finding + fix: `gpu_image_unavailable` → fix the registry `image_uri`
  (image is pulled, never built locally); `lock_stale_vs_yaml` → `workflow/envs/processing.yaml`
  newer than the lock — run `Processing-MuAgent regenerate-locks` (needs `pip install '.[dev]'`)
  and commit; `platform_unsupported` → CPU env is linux-only (use a linux host/container);
  `provision_failed`/`import_failed` → relay the stderr tail. Brand-new machine never
  bootstrapped → `Execution-MuAgent init-machine --processing-repo <repo>` first.

## Supervisor / sentinels

- **Per-stage specs not written** — written automatically by `executor plan-review`. If
  `internal/stage_meta/` is missing/empty, re-run `executor plan-review --config $CFG`.
  Internal state — don't surface unless asked.
- **`hpc-status` shows "Supervisor: not running" with a RUNNING/PENDING scheduler state** —
  the daemon died but the job is alive (stalls won't be auto-cancelled). Restart:
  `executor supervisor-restart --config $CFG` (resumes the watch loop without resubmitting).
  Tell the user what happened.
- **Daemon crashes on a `KillUserProcesses=yes` site** — the SSH session ending killed it;
  the job keeps running but unprotected. Use `supervisor-restart` on reconnect; going
  forward, run `submit` inside `tmux`/`screen`.
- **One-shot `hpc-status` shows "review gate awaiting approval" before any progress** — a
  stale `awaiting_approval` sentinel from a prior run. Fix: (1)
  `squeue -u $USER | grep "pma_head_job_$(basename <run_dir>)"` → `scancel` each; (2)
  `rm internal/proposals/<stage>.awaiting_approval`; (3) `executor submit` and report the
  next one-shot `hpc-status`.
- **Tempted to monitor a long-running job yourself** — don't. The daemon is the sole
  monitor; read one-shot `executor hpc-status`. Re-poll only via a non-blocking scheduled
  wakeup ([`hpc_monitoring.md`](hpc_monitoring.md)) — never a blocking loop or `tail -f | grep`.
