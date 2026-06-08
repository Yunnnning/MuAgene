# QC threshold revision at `post_qc_review`

Use this procedure whenever the user asks to adjust QC filtering or doublet thresholds while the pipeline is paused at the **QC review checkpoint** (`post_qc_review`).

Canonical user-facing report after re-run: `deliverables/checkpoint/qc_review/qc_review_<run_name>.md`. The rendered HTML report is `deliverables/checkpoint/qc_review/qc_summary_<run_name>.html`.

---

## Plan vs live parameters

`executor revise` does **not** update the preprocessing plan.

| Artifact | Updated by `revise`? | Role after QC revision |
|---|---|---|
| `internal/parameters.yaml` | Yes | Live source of truth; stages read overrides here on re-run |
| `internal/artifacts/p2_plan/preprocessing_plan.json` | No | Frozen P2 plan snapshot |
| `deliverables/pre_run/summary/plan_review.md` | No | Plan-review appendix may still show original P2 defaults |
| `deliverables/checkpoint/qc_review/qc_review_<run>.md` | Yes, after re-run + `propose post_qc_review` | User-facing summary of thresholds actually applied |

Do **not** re-run `plan-review` or tell the user the plan file was updated unless they explicitly ask to refresh plan documentation. After a QC revision, the QC review report is the authoritative user-facing summary of applied thresholds.

---

## Allowed side effects beyond CLI

Do not edit `parameters.yaml`, `state.yaml`, biological-context files, or checkpoint sentinels by hand.

Exception for this procedure only: you may delete the stale artifact files listed below under `internal/artifacts/` so Snakemake re-executes affected stages. Do not delete `atac_fragments_cbf_chrnorm.tsv.gz*`; it is expensive to regenerate and is intentionally preserved.

---

## Procedure for HPC runs

Use this when `execution.mode` is `pbs` or `slurm`.

1. **Update parameters.** For each changed value, run:

   ```bash
   executor revise <stage> <stage>.<param>=<value> --config $CFG --rationale "<user's reason>"
   ```

   Valid QC revision stages here are `s1_rna_qc`, `s2_atac_qc`, and `s3_doublets`.

2. **Delete stale artifacts** so Snakemake re-runs the affected stages:

   - **S1 revised:** delete `internal/artifacts/s1_rna_qc/rna_qc.h5ad` and `internal/artifacts/s1_rna_qc/qc_summary.json` (the JSON is the stage-done marker — deleting it causes stage_progress to show S1 as incomplete while it re-runs), plus all S3 artifacts below.
   - **S2 revised:** delete `internal/artifacts/s2_atac_qc/atac_qc.h5ad`, `internal/artifacts/s2_atac_qc/atac_snap.h5ad`, and `internal/artifacts/s2_atac_qc/qc_summary.json`; keep `internal/artifacts/s2_atac_qc/atac_fragments_cbf_chrnorm.tsv.gz*`, plus delete all S3 artifacts below. **Note:** `atac_snap.h5ad` may already be absent if `post_qc_review` was previously approved — the approval cleanup deletes it. The deletion step is safe if the file is missing.
   - **S3 revised:** delete `internal/artifacts/s3_doublets/rna_post_doublet.h5ad`, `internal/artifacts/s3_doublets/atac_post_doublet.h5ad`, `internal/artifacts/s3_doublets/calls.parquet`, `internal/artifacts/s3_doublets/joint_barcodes.txt`, and `internal/artifacts/s3_doublets/overlap_summary.json`.
   - Any S1 or S2 revision invalidates S3. Always delete all five S3 artifacts when S1 or S2 changes.

3. **Approve stages that must re-run.**

   - Approve every revised stage:

     ```bash
     executor approve <stage> --config $CFG
     ```

   - If S1 or S2 was revised, also approve S3 even if you did not run `revise s3_doublets`:

     ```bash
     executor approve s3_doublets --config $CFG
     ```

   - Do **not** pass `--auto-approve` to `submit`. It can refresh sentinels and trigger spurious re-runs of already-complete stages.

4. **Submit and monitor.**

   Cancel stale head jobs first. An old head job can recreate sentinels you just cleared and confuse the new run:

   ```bash
   squeue -u $USER | grep pma_head
   # scancel <JOBID> for each listed pma_head job
   ```

   Then submit and monitor:

   ```bash
   source deliverables/pre_run/config/hpc.env
   executor submit --config $CFG --executor pbs|slurm
   executor hpc-status --watch --config $CFG
   ```

   Use `hpc-status --watch` only. Never substitute `tail -f | grep`.

5. **Regenerate QC reports.** The inferred submit target is an execute rule such as `s3_doublets_execute`; the head job exits after S3 and does not run the local propose rule. After S3 completes successfully, run:

   ```bash
   executor propose post_qc_review --config $CFG
   ```

   This rewrites `qc_review_<run>.md`, `qc_summary_<run>.html`, and checkpoint figures.

6. **Surface the updated report.**

   Read `deliverables/checkpoint/qc_review/qc_review_<run_name>.md` and relay it **verbatim**. Point the user to `qc_summary_<run_name>.html` for the rendered report. Ask whether to approve QC and continue to dimensionality reduction and clustering, or revise again.

Do not surface `internal/proposals/post_qc_review.yaml` in place of the QC review markdown after report regeneration.

---

## Procedure for local runs

Use this when `execution.mode` is `local`.

1. Run the same parameter update step as the HPC procedure.
2. Delete the same stale artifacts as the HPC procedure.
3. Approve every revised stage, and also approve `s3_doublets` when S1 or S2 changed.
4. Re-run execute steps through S3:

   ```bash
   executor run --config $CFG --target s3_doublets_execute
   ```

   Do not use `--auto-approve`.

5. Regenerate QC reports:

   ```bash
   executor propose post_qc_review --config $CFG
   ```

6. Read `deliverables/checkpoint/qc_review/qc_review_<run_name>.md` and relay it **verbatim**. Ask whether to approve QC or revise again.

---

## Branch note

On paired runs, union doublet removal policy is confirmed at `post_qc_review`. It is not a separate S3 user gate. Threshold revisions at S1, S2, or S3 use this procedure only.

---

## Common revise keys

| Stage | Common keys |
|---|---|
| `s1_rna_qc` | `total_counts_min`, `total_counts_max`, `n_genes_min`, `n_genes_max`, `pct_counts_mt_max`, `pct_counts_ribo_max` |
| `s2_atac_qc` | `n_fragments_min`, `n_fragments_max`, `tss_enrichment_min`, `tss_enrichment_max`, `nucleosome_signal_max`, `frip_min` |
| `s3_doublets` | `rna_doublet_score_threshold`, `atac_doublet_score_threshold` |

## Stage-done sentinels

`internal/artifacts/s1_rna_qc/qc_summary.json` and `internal/artifacts/s2_atac_qc/qc_summary.json` are the stage-done markers used by `stage_progress.py` and by Execution-MuAgent for output verification. Both files are written by the respective stage and persist after the `post_qc_review` approval cleanup (which only removes the large h5ads). When revising S1 or S2, delete these JSON files alongside the h5ads so the pipeline correctly reflects that those stages need to re-run.

## QC report cell-count columns

After re-run, `qc_review_<run>.md` and `qc_summary_<run>.html` show:

- **Before / retained / removed** — one summary block per modality at the top of the RNA and ATAC sections (no intermediate ATAC waterfall).
- **cells removed\*** — order-independent exclusive counts from `cells_removed_per_metric` in each stage's `qc_summary.json`. Each per-metric row counts cells failing only that threshold while passing all others in the group. Cells failing multiple thresholds appear under `multiple_metrics`. `total_removed` is the overall removal count for that stage (for ATAC, includes FRiP failures among core-metric passers).
