---
name: downstream_dimred_clustering
domain: downstream
purpose: Document the unattended finish batch — normalization, dimensionality reduction, neighbor graph, clustering, UMAP — that runs after post_qc_review approval up to manifest.
activation: post_qc_review approved; the finish batch (S4–S8 + manifest) is running
inputs: [post_qc_review.approved, internal/parameters.yaml, internal/stage_meta/*.yaml]
outputs: [deliverables/results/<processed h5mu|h5ad>, deliverables/figures/<umap>]
calls_tools: [submit, run, status]
reads_contracts: [stage_meta, latest_snapshot]
writes_state: []
handoff: { next: completion_handoff, when: manifest complete, on_error: troubleshooting }
---

# Downstream — dimensionality reduction & clustering (S4–S8, unattended)

Entered after `post_qc_review` is approved **and** `qc_handoff` has already run (see
[`qc_review_and_revise.md`](qc_review_and_revise.md) — `qc_handoff` runs at the approval
step, not here). This skill is entered when the finish batch is submitted:
`executor submit --config $CFG --executor slurm` (target `all`; Snakemake skips
`qc_handoff` since its outputs already exist) or `executor run` locally. From here to
`manifest` there is **no user checkpoint** — this skill documents what runs so you can
report progress and recognise a healthy finish.

## What runs (per modality, branch-aware)

| Stage | What it does |
|---|---|
| `qc_handoff` | **Already ran at QC approval** — wrote `deliverables/qc/post_qc_<run>.h5mu` (+ manifest) and **deleted** the internal `s3_doublets/{rna,atac}_post_doublet.h5ad` (the h5mu is the canonical post-QC store). Snakemake skips it here. |
| S4 `s4_rna_norm` | RNA normalization (log-normalize, `target_sum=1e4`) + HVG (`seurat_v3` on counts). **Reads RNA from the post-QC h5mu** (rna mod), not the internal h5ad. |
| S5 `s5_atac_spectral` | ATAC spectral embedding via SnapATAC2 (+ flexible peak/feature export). **Reads ATAC from the post-QC h5mu** — rebuilds a snap-native working file from the atac mod's fragments. |
| S6 `s6_neighbors` | PCA (RNA) + neighbor graph — RNA PCA+neighbors; ATAC KNN on the S5 spectral embedding. |
| S7 `s7_clustering` | Leiden clustering at **fixed per-modality resolutions** (values live in `executor/defaults.py` → `s7_clustering`; change them only at plan review via `revise s7_clustering …`). |
| S8 `s8_umap` | UMAP per modality + final write: `processed_<run>.h5mu` (paired) or separate `*_processed.h5ad`. **Hard stop.** |
| `manifest` | Writes `run_manifest.json` and finalizes `deliverables/results/` (localrule). |

For `rna_only`/`atac_only`, the irrelevant per-modality stages are dropped from the plan and
DAG automatically. On `paired`, S7 labels are diagnostic (UMAP only); on
`separate`/single-modality they are the final `leiden_rna`/`leiden_atac` labels.

## Your role here

- **Do not** intervene — no gate, no revise. S7 resolutions are fixed by the plan; a
  different value requires `revise s7_clustering …` back at plan review, not here.
- **HPC:** follow **report-and-repoll** ([`hpc_monitoring.md`](hpc_monitoring.md)); report
  only when the `State:` fingerprint changes. Stop when `monitor.pid` is gone or `manifest`
  is complete.
- **Local:** the finish batch runs straight through under `executor run`.
- On any runtime failure → [`troubleshooting.md`](troubleshooting.md) ("a stage execute
  fails at runtime").

When `manifest` completes → [`completion_handoff.md`](completion_handoff.md).
