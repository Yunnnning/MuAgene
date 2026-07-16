# Example run showcase — 10k PBMC Multiome

Output copy of MuAgene preprocessing deliverables from a paired multiome
run on the public
[10k Human PBMCs, Multiome v1.0, Chromium X](https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-x-1-standard-2-0-0)
dataset.

## What is included

| Path | Contents |
|---|---|
| `deliverables/plan/` | Plan review markdown/HTML, biological context, sanitized `run.yaml` |
| `deliverables/qc/` | QC review markdown/HTML and `post_qc_manifest.json` |
| `deliverables/figures/` | data exploration, QC filtering, and UMAP figures (PNG/PDF) |
| `deliverables/results/` | `run_manifest.json`, `layout.json`, review notebook/script |

Final cell count after QC and doublet removal: **7,531** jointly retained cells.

Referenced but not shipped:

- `deliverables/results/processed_pbmc10k_multiome.h5mu`
- `deliverables/qc/post_qc_pbmc10k_multiome.h5mu`

## Quick look

- Plan: [`deliverables/plan/plan_review_pbmc10k_multiome.md`](deliverables/plan/plan_review_pbmc10k_multiome.md) /
  [`plan_summary_pbmc10k_multiome.html`](deliverables/plan/plan_summary_pbmc10k_multiome.html)
- QC: [`deliverables/qc/qc_review_pbmc10k_multiome.md`](deliverables/qc/qc_review_pbmc10k_multiome.md) /
  [`qc_summary_pbmc10k_multiome.html`](deliverables/qc/qc_summary_pbmc10k_multiome.html)
- UMAPs: [`deliverables/figures/s8_umap_rna_by_leiden.png`](deliverables/figures/s8_umap_rna_by_leiden.png),
  [`deliverables/figures/s8_umap_atac_by_leiden.png`](deliverables/figures/s8_umap_atac_by_leiden.png)
