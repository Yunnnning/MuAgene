# MuAgene state model

Every persistent state file MuAgene reads or writes, with its owner (writer),
consumer (reader), and lifecycle. This is the catalog the prompts/skills reference
instead of describing state inline. Paths are relative to a run directory unless
marked machine-level.

## Golden rule
All run state is mutated **only** through the `executor` CLI (Processing) or the
`Execution-MuAgent` CLI. Prompts never hand-edit these files. Biological context uses
its deterministic mapper/writer API; QC revision and cleanup remain CLI-owned.

## Run state — written by Processing-MuAgent (`executor`)

| Path | Writer (command) | Reader | Lifecycle |
|------|------------------|--------|-----------|
| `internal/parameters.yaml` | `init`, `declare-branch`, `configure-execution`, `revise`, stage runs (echo-back) | every stage; `plan_assembler.overlay_plan`; reports | Live source of truth for parameters; provenance per key `{value, source, confidence, rationale, method, revision_of}`. |
| `internal/artifacts/p2_plan/preprocessing_plan.json` | `plan-review` (assemble_plan + write_plan) | all stages (the default layer; parameters.yaml overrides it) | Frozen plan snapshot; re-rendered if the plan is rebuilt. |
| `internal/proposals/<stage>.yaml`, `<stage>.awaiting_approval` | Snakemake `*_propose` localrules | the agent (relayed to the user at a gate) | Created when a gate is reached; consumed at approval. |
| `internal/checkpoints/<gate>.approved` | `approve <gate>` | Snakemake execute rules (hard gate) | Written on approval; presence unblocks downstream. Gates: `plan_review`, `post_qc_review`. |
| `internal/artifacts/<stage>/*` | the stage `run()` | downstream stages; reports | Durable markers preserve DAG/status state. Regenerable QC files are cleaned on QC approval, S3 post-doublet files after a successful handoff, and S4–S8 working files by `finish-cleanup`. Authoritative deletion sets live in `executor/cleanup.py` and `executor/stages/qc_handoff.py`; prepared ATAC fragments are retained by default for integration. |
| `internal/stage_meta/<stage>.yaml`, `internal/stage_meta/head_job.yaml` | `plan-review` / `submit` (specs) | **Execution-MuAgent** (the spec contract) | Per-stage science intent + resources + I/O + progress_timeout_hint; head-job submission spec. |
| `internal/log.jsonl` | `log_event` (all stages) | debugging / audit | Append-only event log. |
| `deliverables/plan/config/run.yaml` | `init` (canonical copy) | every CLI call (`--config`) | The run config; canonical path after `init`. |
| `deliverables/plan/config/biological_context.md` | `context_mapper.write_report` | P1 context stage | Blank template or filled; the only non-CLI write path (still deterministic). |
| `deliverables/plan/config/site.config` | `configure-execution` (HPC) | **Execution-MuAgent** | Platform description (scheduler, partition/account, device+GPU routing, env identity, `environments:` recipe). |
| `deliverables/plan/config/hpc.env` | `configure-execution` (HPC) | sourced before `submit` | Shell exports for the cluster. |
| `deliverables/plan/plan_review_<run>.md`, `plan_summary_<run>.html` | `plan-review` | user (checkpoint #1) | Re-rendered deterministically from canonical sources. |
| `deliverables/qc/qc_review_<run>.md`, `qc_summary_<run>.html` | `post_qc_review` propose | user (checkpoint #2) | QC summary of filters actually applied. |
| `deliverables/qc/post_qc_<run>.h5mu`, `peaks_<run>.bed`, `post_qc_manifest.json` | `qc_handoff` | downstream consumer | Post-QC handoff bundle — schema `muagene.post_qc_handoff/1` (see `post_qc_manifest.schema.json`). |
| `deliverables/results/processed_<run>.h5mu` / `*_processed.h5ad` | S8 | user / downstream | Final per-modality processed output. |
| `deliverables/results/run_manifest.json` | `manifest` | user / downstream | Preprocessing handoff manifest (v1.0.1). |

## Shared execution state (under `internal/hpc_monitor/`)

| Path | Writer | Reader | Lifecycle |
|------|--------|--------|-----------|
| `latest_snapshot.json` | `execute-spec` preflight and monitor watch loop | **Processing-MuAgent** (`hpc-status`) | Structured run contract. Pre-submit failures write `PRE_SUBMIT_FAILED`; active monitoring refreshes scheduler/log state, findings, actions, and recheck cadence. |
| `latest_report.md` | monitor | explicit operator `report` command | Debug/audit only; never parsed by Processing or used for normal user status. |
| `monitor.pid` | Processing `submit` / `supervisor-restart` | Processing (liveness signal); Execution removes | Written when Processing starts the supervisor; removed when Execution monitoring exits or Processing terminates/replaces the supervisor. |
| `latest_submission.json` | `execute-spec` | `resume-monitor`; Processing status | Current submission context used to resume monitoring without resubmitting. |
| `submissions.jsonl`, `execution_manifest.jsonl` | `execute-spec` / `resume-monitor` | audit | Append-only submission/registration logs. |
| `scripts/<stage>_<ts>.sh` | submit (rendered) | scheduler (sbatch) | Rendered submission script per submit. |

## Machine-level state — written by Execution-MuAgent (under `~/.muagene/`)

| Path | Writer | Reader | Lifecycle |
|------|--------|--------|-----------|
| `machine.config` | `init-machine` | provision/validate/execute-spec | Per-host facts (manager, container runtime, GPU image + pinned image_uri, policy, env names, processing-repo path). |
| `env_state.json` | `provision-env` (record) | execute-spec preflight, `validate-env` | Per-device env fingerprints (`lock:<env>` -> lock sha256; `container:<sif>` -> image_uri). `missing`/`stale` -> auto-provision (policy=auto) or fail loud (manual). |
| `images/muagene-gpu.sif` | `provision-env --device gpu` (pull-only) | GPU child jobs | Pinned GPU container; pulled, never built locally. |

## Monitor states (not finding codes)
`HEALTHY`/`RECOVERED` -> `SUSPECT` -> `INVESTIGATING` -> `CONFIRMED_DEAD` | `FS_HANG` | `RECOVERED`; plus `DONE` (clean finish) and `KILLED`. `confirmed_dead_reason` rides in `kill_action`. The actionable signals Processing consumes are the `findings` (see `findings.yaml`), not the raw state.
