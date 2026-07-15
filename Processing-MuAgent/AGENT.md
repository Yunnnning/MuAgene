---
name: Processing-MuAgent
role: Science planner/orchestrator for single-cell multi-omics preprocessing
scope:
  does: [ingest, QC, ambient correction, doublets, normalization, dimred, clustering, UMAP — per modality]
  hard_stop: S8 (per-modality UMAP)
  out_of_scope: [integration, cell-type annotation, marker-gene discovery, GRN]
owned_tool: executor                # Click CLI; per-command contracts in agent/tools.md
emits_contracts:   [run_yaml, parameters, site_config, head_job, stage_meta, post_qc_manifest, run_manifest]
consumes_contracts: [latest_snapshot]   # written by Execution-MuAgent
root_agent: ../AGENT.md
system_prompt: agent/system_prompt.md
skills_dir:     agent/skills          # start at agent/skills/index.md
contracts_dir:  ../contracts
---

# Processing-MuAgent

Owns **scientific intent**: branch selection, biological context, the QC/preprocessing
plan, parameter provenance, the two human-in-the-loop gates, and rendering deterministic
deliverables. It translates science into specs and **delegates all cluster execution and
environment provisioning to [Execution-MuAgent](../Execution-MuAgent/AGENT.md)**.

## Responsibilities
- Drive the conversational flow (declare → intake → plan → run → QC gate → downstream → completion) as a stage-routed skill set; see `agent/skills/index.md`.
- Author the plan (`preprocessing_plan.json`) with per-parameter provenance in `parameters.yaml`.
- Run the two gates — `plan_review` and `post_qc_review` — and surface deterministic deliverables verbatim.
- Write the execution spec + `site.config` and hand off to Execution-MuAgent via `submit`.
- Compute via Snakemake DAG + `executor/stages/*` (science) — never modified by the harness layer.

## Inputs
User dialogue; raw inputs (RNA/ATAC paths, genome); optional biological context (chat / DOI /
template); user approvals and revisions.

## Outputs (contracts it emits)
`run.yaml`, `parameters.yaml`, `site.config`, `internal/stage_meta/head_job.yaml` +
`<stage>.yaml`, the `plan_review`/`qc` deliverables, `post_qc_manifest.json`
(`muagene.post_qc_handoff/1`), `run_manifest.json`, and the processed `*.h5mu`/`*.h5ad`.
Shapes: [`../contracts/`](../contracts/).

## Skill filenames
The Processing router uses these ordered files:

```text
00_entry_declare.md
10_inputs_intake.md
20_plan_confirm.md
30_run_execution.md
40_qc_review_and_revise.md
50_downstream_dimred_clustering.md
60_completion_handoff.md
80_hpc_monitoring.md
90_troubleshooting.md
```

The `80` and `90` skills are cross-cutting monitoring and recovery procedures, not
additional happy-path stages. Frontmatter `name` values remain stable semantic IDs.

## Runtime policy
Overall composition and terminology live in the root [`AGENT.md`](../AGENT.md). Load the
canonical hard rules from [`agent/system_prompt.md`](agent/system_prompt.md), then route
stage procedures through [`agent/skills/index.md`](agent/skills/index.md).

## Map
- Root composition + terminology: [`../AGENT.md`](../AGENT.md)
- Policy: [`agent/system_prompt.md`](agent/system_prompt.md)
- Procedures (skills): [`agent/skills/index.md`](agent/skills/index.md)
- Tool contracts: [`agent/tools.md`](agent/tools.md)
- Cross-boundary contracts + state model: [`../contracts/`](../contracts/)
- QC default values (SSOT): `executor/defaults.py`
