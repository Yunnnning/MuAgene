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
hard_rules: [no-invented-values, confirm-exec-mode-once, state-via-cli-only, verbatim-output, no-silent-overrides, stop-at-S8]
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
- Drive the conversational flow (declare → intake → plan → run-with-checkpoints); see `agent/skills/`.
- Author the plan (`preprocessing_plan.json`) with per-parameter provenance in `parameters.yaml`.
- Run the two gates — `plan_review` and `post_qc_review` — and surface deterministic deliverables verbatim.
- Write the execution spec + `site.config` and hand off to Execution-MuAgent via `submit`.
- Compute via Snakemake DAG + `executor/stages/*` (science) — never modified by the harness layer.

## Inputs
User dialogue; raw inputs (RNA/ATAC paths, genome); optional biological context (chat / DOI /
template); user approvals and revisions.

## Outputs (contracts it emits)
`run.yaml`, `parameters.yaml`, `site.config`, `internal/stage_meta/head_job.yaml` +
`<stage>.yaml`, the `plan_review`/`qc_review` deliverables, `post_qc_manifest.json`
(`muagene.post_qc_handoff/1`), `run_manifest.json`, and the processed `*.h5mu`/`*.h5ad`.
Shapes: [`../contracts/`](../contracts/).

## Constraints
Never invent paths/values/genes; record state only via the `executor` CLI; confirm execution
mode once (`execution.user_confirmed`); never submit or monitor cluster jobs itself; surface
executor output verbatim; **hard-stop at S8**.

## Failure modes
Missing input → ask/wait. Pairing ambiguous → relay the raised error (no silent retry).
Stale CPU lock → fail loud (`lock_stale_vs_yaml`). Destructive `revise` at `post_qc_review`
→ run the binding-constraint diagnosis + `--dry-run` and confirm first
([`agent/skills/qc_review_and_revise.md`](agent/skills/qc_review_and_revise.md)).

## Map
- Policy + entry point: [`agent/system_prompt.md`](agent/system_prompt.md)
- Procedures (skills): [`agent/skills/index.md`](agent/skills/index.md)
- Tool contracts: [`agent/tools.md`](agent/tools.md)
- Cross-boundary contracts + state model: [`../contracts/`](../contracts/)
- QC default values (SSOT): `executor/defaults.py`
