---
name: plan_confirm
domain: plan
purpose: Kick off planning compute, render the plan with an intro, resolve the marker-gene and QC-threshold questions, then run the plan_review gate (#1).
activation: run scaffolded (init done) + branch declared + exec-mode confirmed; plan_review not yet approved
inputs: [deliverables/plan/config/run.yaml, deliverables/plan/config/biological_context.md, internal/artifacts/s0_ingest/validation_report.json]
outputs: [deliverables/plan/plan_review_<run>.md, deliverables/plan/plan_summary_<run>.html, internal/stage_meta/*.yaml, plan_review.approved]
calls_tools: [run, submit, hpc-status, plan-review, revise, approve]
reads_contracts: [stage_meta, parameters]
writes_state: [parameters.yaml, plan_review.approved]
handoff: { next: run_execution, when: plan_review approved, on_error: troubleshooting }
---

# Plan confirmation — the plan_review gate (#1)

Entered after [`10_inputs_intake.md`](10_inputs_intake.md) has scaffolded the run, declared the
branch, and confirmed execution mode. This skill runs planning compute, presents the plan,
and drives the first human gate. **Nothing heavy runs until the user approves.**

## 1. Kick off planning (P1 → S0 → gate-arming)

The planning target is **`plan_review_propose`** (auto-inferred when `--target` is omitted).
It depends on `s0_ingest_execute`, so one invocation runs P1 → S0 → plan assembly →
gate-arming and arms `plan_review` at the end. Do **not** target `s0_ingest_execute` alone —
that stops one rule early and leaves the gate unarmed.

- **HPC** (`execution.mode = slurm`):
  ```
  source deliverables/plan/config/hpc.env
  executor submit --config $CFG --executor slurm     # auto-infers plan_review_propose; S0 needs 100+ GB → compute node
  executor hpc-status --config $CFG                   # one-shot, then report-and-repoll
  ```
  Then follow [`80_hpc_monitoring.md`](80_hpc_monitoring.md) until `plan_review` is
  `awaiting_approval`.
- **Local** (`execution.mode = local`):
  ```
  executor run --config $CFG --target plan_review_propose       # ~30s on small inputs
  ```

Do **not** retry logic errors (pairing ambiguous, declared-vs-detected mismatch, missing
index, S0 OOM) — see [`90_troubleshooting.md`](90_troubleshooting.md). A declared `paired`
run that cannot be validated stops in S0; obtain an explicit branch/input resolution
before rerunning.

## 2. Render the plan with an intro paragraph

**Prerequisite:** `executor plan-review` is a *renderer* that requires the planning
compute (`plan_review_propose`) to have finished and produced
`internal/artifacts/p2_plan/preprocessing_plan.json` and
`internal/artifacts/s0_ingest/validation_report.json`. Do **not** run it before the
Snakemake planning job has completed — the CLI will now refuse and the command would
otherwise emit placeholder deliverables and a false `plan_review.awaiting_approval`
sentinel.

1. `executor plan-review --intro-context --config $CFG` — prints JSON (sample metadata,
   cell counts, barcode matching). Write nothing yet.
2. Write a 100–150-word intro paragraph from that data. Cover organism, tissue,
   platform/assay, the aim (QC → doublet removal → dimred → clustering), raw cell counts per
   modality, and the barcode-matching result. Smooth prose, no bullet points, no stage codes
   or internal filenames, no rounded/omitted numbers.
   - **Paired-candidate compatibility check (only when `workflow_branch = paired`):** if the
     barcode check found no direct/subset match (`pairing_confidence` not "high" or
     `pairing_status` not "paired"), diagnose rather than just flag: read `pairing_ladder`
     (which rungs were attempted + why each failed), cross-check `rna_filtered_status` /
     `atac_barcodes_source` and `rna_raw_n_barcodes` vs `rna_n_cells` (a raw ATAC file
     explains near-zero overlap), and read `run.yaml` for mismatched sources. Put the most
     likely root cause + concrete fix (e.g. use the filtered matrix, supply a
     `barcode_translation_path`) into the intro. Skip for `rna_only`/`atac_only`/`unpaired`.
3. `executor plan-review --intro "<paragraph>" --config $CFG` — re-renders BOTH
   `plan_review_<run>.md` and `plan_summary_<run>.html` with the intro prepended, persists
   the intro (pass `--intro` once), and writes the per-stage specs to `internal/stage_meta/`.

## 3. Marker-gene check — mandatory question when ambient correction is planned

If the plan keeps ambient RNA correction (`s1a_ambient.method != none`) and the rendered
"Marker gene expression check" item is still `not set`, you **must** ask before any approval
(escalate to *strongly recommended* when `qc_explore` median rho is high):

> The plan runs ambient RNA correction. I recommend checking marker-gene expression
> *before vs after* correction. Please give me 5–10 marker gene symbols to visualise, or
> tell me to **defer** this to QC review, or to **skip** it.

**Never invent, suggest, or look up gene names** — canonical rule in
[`40_qc_review_and_revise.md`](40_qc_review_and_revise.md). Record the user's one explicit choice:
- **provide genes** → `executor revise s1a_ambient "marker_genes=[g1, g2, ...]" --config $CFG --rationale "Marker genes provided at plan review"` (plotted automatically during S1a).
- **defer** → carry `--defer-marker-genes` on the approve call (`--marker-genes defer` on `submit`).
- **decline** → carry `--skip-marker-genes` on the approve call (`--marker-genes skip` on `submit`).

If `s1a_ambient.method == none`, skip this question. The executor refuses to approve while
the decision is unresolved.

## 4. QC threshold confirmation — mandatory, after the marker-gene step

This is the **plan_review home for QC threshold revision** — the user can set, adjust, pin, or
skip any QC threshold here, *before* any stage runs. The "QC strategy" item shows
`[? needs confirmation]`. Ask:

> The default MAD-based thresholds are in the plan appendix histograms. Keep the defaults,
> adjust one or more, pin an exact value (e.g. RNA `n_genes` lower bound = 300), or skip a
> metric entirely?

- **Keep defaults** → go to step 5.
- **Adjust / pin / skip** → `executor revise s1_rna_qc <param>=<value> --config $CFG --rationale "<reason>"` (or `s2_atac_qc`). Revise the input knobs, or pin a bound with its `*_override` key — the common keys, the `*_override` pinning table, and the skip recipes are the shared reference in [`40_qc_review_and_revise.md`](40_qc_review_and_revise.md) (they apply at both gates). The binding-constraint diagnosis there (a MAD bound is `max(MAD, floor)` / `min(MAD, ceiling)`, so a non-binding knob has no effect) is just as useful here for picking the right knob.

**`revise` at plan_review is non-destructive.** Nothing has run yet, so it only updates
`internal/parameters.yaml` and re-renders the plan deliverables + the S0 QC-exploration
preview (recomputed projected removals) — **no artifact deletion and no stage re-run** (unlike
post_qc_review, where `revise` deletes downstream QC artifacts; see
[`40_qc_review_and_revise.md`](40_qc_review_and_revise.md)). So no dry-run guardrail is needed here:
issue the `revise`, re-surface the regenerated plan, and re-ask. The approve / revise / abort
loop is described once in step 5.

## 5. Approve / revise / abort

- **Approve** → `executor approve plan_review --config $CFG --note "approved after review"`
  (+ `--defer-marker-genes`/`--skip-marker-genes` to match step 3 when no genes were given;
  on HPC the same is `--marker-genes defer|skip` on `submit --auto-approve`).
- **Revise** → `executor revise <stage> <param>=<value> --config $CFG --rationale "<reason>"`
  (stage prefix auto-added). While `plan_review` is unapproved, `revise` auto-regenerates the
  plan deliverables; the stage returns to `awaiting_approval`. Re-surface and re-ask.
- **Abort** → stop; the run dir is intact and resumable on the same config.

## What to surface back

- The **Summary** section of `plan_review_<run>.md`, **verbatim** (the appendix is reference
  detail). Don't paraphrase values. Point the user at `plan_summary_<run>.html`.
- If marker genes were stored, confirm the gene list in one line.
- On approval, hand off to [`30_run_execution.md`](30_run_execution.md) for the QC batch.
