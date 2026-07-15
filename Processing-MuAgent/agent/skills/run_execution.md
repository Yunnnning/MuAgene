---
name: run_execution
domain: execute
purpose: Dispatch local run / HPC submit of the next stage batch, drive the checkpoint loop, and arm the next gate. The execution engine for both the QC batch and the finish batch.
activation: a gate was just approved (plan_review or post_qc_review) and the next batch must run; or a batch is mid-flight
inputs: [deliverables/plan/config/run.yaml, deliverables/plan/config/hpc.env, deliverables/plan/config/site.config, internal/stage_meta/*.yaml]
outputs: [internal/stage_meta/head_job.yaml, internal/proposals/<stage>.yaml, internal/proposals/<stage>.awaiting_approval, internal/hpc_monitor/*]
calls_tools: [run, submit, status, approve]
reads_contracts: [head_job, stage_meta, site_config, latest_snapshot]
writes_state: [head_job.yaml, proposals/<stage>.awaiting_approval]
handoff: { next: qc_review_and_revise OR downstream_dimred_clustering, when: a gate is awaiting_approval or the finish batch is running, on_error: troubleshooting }
---

# Run execution — the dispatch + checkpoint engine

Entered after a gate approval — from [`plan_confirm.md`](plan_confirm.md) for the **QC
batch**, or from [`qc_review_and_revise.md`](qc_review_and_revise.md) for the **finish
batch**. This skill only *moves compute and arms gates*; the science at each gate is owned by
that gate's own skill.

## Execution boundary (hard)

`executor run` is **local-only**; `executor submit` is **cluster-only**. Processing-MuAgent
**never submits or monitors cluster jobs itself** — `submit` writes the head-job spec and
delegates all cluster execution **and** env provisioning to Execution-MuAgent (consuming the
`site.config` written at intake). There is no `run --executor slurm`. Read `execution.mode`
from `parameters.yaml` (ask if missing) to pick the path.

## `submit` mechanics (HPC)

`submit` writes `internal/stage_meta/head_job.yaml`, starts `Execution-MuAgent execute-spec`
as a **background supervision daemon** (submits the head-job, records the job ID to
`execution_manifest.jsonl`, runs the watch loop for the job's lifetime), returns within
~90 s once the job ID is confirmed, and **fails loudly if Execution-MuAgent is absent**.
Processing writes `monitor.pid` when it starts the daemon; Execution removes it when monitoring
exits. (On `KillUserProcesses=yes` sites, run `submit` inside `tmux`/`screen`.)

## HPC run phases + batch staging

| Phase | Stages | Executes on | Owner skill |
|---|---|---|---|
| Context | P1 | login node (localrule) | — |
| S0 ingest (+plan) | S0 | head-job (HPC) / login (local) | [`plan_confirm.md`](plan_confirm.md) |
| Gate #1 | plan_review | login node | [`plan_confirm.md`](plan_confirm.md) |
| QC | S1a → S3 | head-job | this skill |
| Gate #2 | post_qc_review | — | [`qc_review_and_revise.md`](qc_review_and_revise.md) |
| Integration handoff | qc_handoff | SLURM cluster job at QC approval (`submit --target qc_handoff`) | [`qc_review_and_revise.md`](qc_review_and_revise.md) |
| Finish | S4 → S5 → S6 → S7 → S8 → manifest | head-job | [`downstream_dimred_clustering.md`](downstream_dimred_clustering.md) |

After plan-review approval, `source deliverables/plan/config/hpc.env`, then:
- **QC batch:** `executor submit --config $CFG --executor slurm --auto-approve --auto-approve-except post_qc_review`
- **After QC approval:** the agent immediately submits the separate `qc_handoff` target (see [`qc_review_and_revise.md`](qc_review_and_revise.md)). After it verifies the handoff and the user confirms the finish batch, `executor submit --config $CFG --executor slurm` runs S4→S8→manifest (target `all`; Snakemake skips the completed handoff).

Each gated phase's head-job target is the **gate-arming `*_propose` localrule** (e.g.
`post_qc_review_propose`), not the phase's last execute stage. Snakemake pulls every execute
stage in the phase as a dependency and runs the propose localrule last, so one submission
runs the whole phase **and** arms the gate (`<stage>` → `awaiting_approval`). The finish
phase has no gate → target `all`. You never run `propose` by hand to surface a gate.

## The checkpoint loop

1. `executor status --config $CFG` — find the stage that is `awaiting_approval`.
2. Read `internal/proposals/<stage>.yaml` (and any linked summary, e.g.
   `deliverables/qc/qc_review_<run>.md`).
3. Route by stage via [`index.md`](index.md): `post_qc_review` →
   [`qc_review_and_revise.md`](qc_review_and_revise.md); other gated stages → approve, or
   `executor revise <stage> <key>=<value> --config $CFG` then re-surface.
4. Approve → `executor approve <stage> --config $CFG`; re-submit/re-run the next batch.
5. Continue until `manifest` completes → [`completion_handoff.md`](completion_handoff.md).

The only human checkpoints are `plan_review` and `post_qc_review`. The handoff-to-finish
confirmation is conversational agent policy, not another Snakemake gate. `s7_clustering`
uses fixed planned resolutions and is not a checkpoint.

## Monitoring + signals

After any `submit`, follow **report-and-repoll** — canonical procedure in
[`hpc_monitoring.md`](hpc_monitoring.md) (never a blocking loop or `tail -f`).
`monitor.pid` removal means the daemon stopped (this phase's compute is over); treat it and
the gate sentinel as **independent** signals — drive the next checkpoint on **either**.

## Job naming

Cluster jobs are `pma_{stage}_{run_name}` (head) / `pma_{rule}_{run_name}` (child), where
`run_name = basename(run_dir)`. **Always filter squeue by run name** — never bare
`grep pma_head` (matches every concurrent sample):
`squeue -u $USER | grep "pma_head_job_$(basename <run_dir>)"`.

## Local mode

`executor run --config $CFG` (no `--auto-approve`) runs every stage whose `.approved`
sentinel exists, stopping at the first missing one. Drive the same checkpoint loop; re-run
after each approval. Local QC-revision re-run details: [`qc_review_and_revise.md`](qc_review_and_revise.md).
