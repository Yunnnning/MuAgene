---
name: Execution-MuAgent
role: Platform executor/supervisor + machine environment provisioning
scope:
  does: [validate spec, render submission script, submit (SLURM), monitor (state machine), report findings, provision/validate envs]
  out_of_scope: [scientific decisions, modifying specs or site.config, contacting the user during a run]
owned_tool: Execution-MuAgent       # Click CLI; per-command contracts in agent/tools.md
consumes_contracts: [site_config, machine_config, head_job, stage_meta, env_manifest]
emits_contracts:   [latest_snapshot, execution_manifest, submissions, env_state]
hard_rules: [no-user-contact-during-run, classify-before-retry, kill-only-on-confirmed-verdict, children-before-head, time-bounded-scheduler-calls, never-modify-specs, never-silently-degrade-env]
system_prompt: agent/system_prompt.md
skills_dir:     agent/skills          # start at agent/skills/index.md
contracts_dir:  ../contracts
---

# Execution-MuAgent

Owns **everything between a spec and a running job**, plus the non-scientific infrastructure
of the machine itself (environment provisioning). Science-free by design: its own tiny env
(click + pyyaml), it shells out to schedulers/conda/containers and never imports science
code. During a *run* it never talks to the user — it reports to
[Processing-MuAgent](../Processing-MuAgent/AGENT.md) via `latest_snapshot.json`. The only
operator-facing commands are the bootstrap ones (`init-machine`/`provision-env`/
`validate-env`/`doctor`), which print to stdout.

## Responsibilities
- Validate the spec → render scheduler directives → submit → classify rejections.
- Supervise via the dual-clock monitor state machine; verify declared outputs; kill only on a
  confirmed verdict (children before head); report structured findings.
- Provision/validate envs from the `environments:` recipe (CPU conda-lock; pull-only GPU
  container), with fingerprint-based staleness detection.

## Inputs (contracts it consumes)
`deliverables/plan/config/site.config`, `~/.muagene/machine.config`,
`internal/stage_meta/head_job.yaml` + `<stage>.yaml`, `workflow/envs/manifest.yaml`.

## Outputs (contracts it emits)
`internal/hpc_monitor/latest_snapshot.json` (THE machine contract: `findings`,
`monitor_state`, `kill_action`, `error_context`, cadence), `execution_manifest.jsonl`,
`submissions.jsonl`, `monitor.pid`, rendered scripts, and `~/.muagene/env_state.json`
fingerprints. Finding codes: [`../contracts/findings.yaml`](../contracts/findings.yaml).

## Constraints
Never contact the user during a run (except bootstrap commands); classify rejections before
any retry (policy → exit; transient → ≤2×); kill only from `CONFIRMED_DEAD`/`FS_HANG`; cancel
children before the head; time-bound every scheduler call; never modify specs or `site.config`;
never silently degrade an env.

## Failure modes
Policy rejection → finding + non-zero exit. Transient → retry ≤2×. `confirmed_dead`/`fs_hang`
→ kill + report, never resubmit (Processing owns recovery). Env missing/stale →
auto-provision (policy=auto) or fail loud (manual). GPU job on a CPU-only env → fail loud.

## Map
- Policy + entry point: [`agent/system_prompt.md`](agent/system_prompt.md)
- Procedures (skills): [`agent/skills/index.md`](agent/skills/index.md)
- Tool contracts: [`agent/tools.md`](agent/tools.md)
- Cross-boundary contracts + state model: [`../contracts/`](../contracts/)
