# processing-MuAgent

Multiome (scRNA-seq + scATAC-seq) preprocessing subagent. Takes filtered or raw 10x Genomics multiome outputs and performs QC, dimensionality reduction, clustering, UMAP-visualisation **per modality**, then **stops** before integration.

## Workflow

Pre-preprocessing phases (run before any QC):

- **P1 Context Extraction** — Biological Context Report (organism, tissue, assay, DOIs)
  + DOI-based prior-analysis extraction.
- **P2 Preprocessing Plan Generation** — holistic `preprocessing_plan.json` for every
  downstream stage; needs approval before execution.

Preprocessing stages:

- **S0 Ingest** — Accepts both Cell Ranger **filtered** and **raw** matrices. Format 
  autodetect (10x h5 / MEX / h5ad / custom), fragments validation (+ tbi), pairing detection 
  (paired vs separate workflow branches), metadata handling (minimal reconstruction when absent). 
- **S1a Ambient RNA correction** — DecontX (filtered counts only) or SoupX
  (raw + filtered) auto-dispatched from `s0` outputs. Pass-through when R / Bioconductor 
  is unavailable.
- **S1 RNA QC** — MAD-derived thresholds on total_counts / n_genes / pct_counts_mt
  + `pct_counts_ribo` ceiling, computed on the decontaminated counts from S1a.
- **S2 ATAC QC** — TSS enrichment + per-cell nucleosome signal (Signac-style
  `mono / nucleosome_free`) + fragment-count MAD via SnapATAC2.
- **S3 Doublets** — Scrublet (RNA, sparse-CSR input, adaptive
  `expected_doublet_rate ≈ 0.0008 × n_cells`) + SnapATAC2 scrublet (ATAC);
  four-way overlap summarised and goal-based removal-policy recommendation
  (union / intersection). Raw calls preserved.
- **S4 RNA norm + HVG** — log-normalize (target_sum=1e4) + HVG (`seurat_v3` on counts).
- **S5 ATAC TF-IDF + LSI and Peak Matrix Export** — Performs TF-IDF normalization and spectral embedding (LSI) on the SnapATAC2 tile matrix (`bin_size=500`, unified with S3). In parallel, the stage attempts to export a feature (cell-by-feature) matrix for ATAC using the following approach:
  1. **Peak Matrix:** If available, loads the peak matrix directly from Cell Ranger ARC (10x) multiome `.h5` files ("arc_h5" mode). If ARC peak matrix is unavailable, generates peaks from fragments using SnapATAC2’s MACS3 integration, then constructs a peak-by-cell matrix ("macs3_from_fragments" mode).
  2. **Fallback:** If both peak matrix approaches fail, exports the verified SnapATAC2 tile matrix ("tile_matrix_fallback" mode).
  The resulting exported matrix is used for downstream integration; the LSI embedding is always preserved for ATAC clustering.
- **S6 Dim reduction + neighbors** — For RNA: applies `sc.pp.scale` (optional), then PCA; the number of principal components (`n_pcs`) is determined using a chord-distance ("elbow") heuristic on the explained-variance curve, capped at `rna_n_pcs_max`. Nearest-neighbors are then computed on the PCA space.
  For ATAC: cells are embedded using SnapATAC2’s LSI (from previous stage), and neighbor graphs are computed directly on the LSI representation.
- **S7 Clustering** — Leiden resolution sweep with per-modality grid, stable-region knee picker; RNA tilt=higher, ATAC tilt=lower.
- **S8 UMAP** — per-modality UMAP; paired → `processed.h5mu`, separate → two `.h5ad`.
- **manifest** — `run_manifest.json` (handoff contract v1.0.0).

## Execution engine

Snakemake: each stage is a `<stage>_propose` + `<stage>_execute` rule pair, with
`<stage>_execute` gated on `checkpoints/<stage>.approved`.

## CLI

**Step 1 — Install** (must run from inside `processing-MuAgent/`):

```bash
cd /path/to/processing-MuAgent
pip install -e .
```

**Step 2 — Edit the example config** at `processing-MuAgent/config/run.example.yaml`.
At minimum set:

```yaml
run_dir:               /path/to/your/output/run_01   # where all outputs will be created
rna_path:              /path/to/filtered_feature_bc_matrix.h5
atac_fragments_path:   /path/to/atac_fragments.tsv.gz
genome_assembly:       GRCh38   # or mm10
```

**Step 3 — Scaffold the run directory** (`init` creates it and copies your config inside):

```bash
# Run from inside processing-MuAgent/ so the --config path resolves correctly.
processing-muagent init --config config/run.example.yaml
```

`init` creates `<run_dir>/` with all internal scaffolding and places two files
for you to review:

- `<run_dir>/deliverables/pre_run/config/run.yaml` — your working config copy
- `<run_dir>/deliverables/pre_run/config/biological_context.md` — organism / tissue / assay template (optional but improves QC threshold selection)

**Step 4 — Run** (point `--config` at the copy that `init` placed inside the run dir):

```bash
# Fully automatic:
processing-muagent run \
    --config <run_dir>/deliverables/pre_run/config/run.yaml \
    --auto-approve

# Or check pipeline status at any point:
processing-muagent status \
    --config <run_dir>/deliverables/pre_run/config/run.yaml
```

**Interactive / checkpoint-by-checkpoint mode:**

```bash
CONFIG=<run_dir>/deliverables/pre_run/config/run.yaml

processing-muagent propose p1_context --config $CONFIG
# review: <run_dir>/internal/proposals/p1_context.yaml
processing-muagent approve p1_context --config $CONFIG

processing-muagent propose p2_plan --config $CONFIG
# review: <run_dir>/deliverables/pre_run/summary/plan_summary.md
processing-muagent approve p2_plan --config $CONFIG

processing-muagent propose plan_review --config $CONFIG
# review: <run_dir>/deliverables/pre_run/summary/plan_review.md
processing-muagent approve plan_review --config $CONFIG

# Then S0 → S1a → S1 → S2 → S3 → S4 → S5 → S6 → S7 → S8:
for STAGE in s0_ingest s1a_ambient s1_rna_qc s2_atac_qc \
             s3_doublets s4_rna_norm s5_atac_lsi \
             s6_dimred s7_clustering s8_umap; do
    processing-muagent propose $STAGE --config $CONFIG
    # review: <run_dir>/internal/proposals/$STAGE.yaml
    processing-muagent approve $STAGE --config $CONFIG
done
```

## Repository layout

```
processing-MuAgent/
├── config/              # example run configurations
├── executor/            # Python implementation (stages, methods, CLI, helpers)
│   ├── stages/          # per-stage scripts S0..S8
│   └── methods/         # named-method helpers (MAD thresholds, resolution sweep, doublet policy)
├── workflow/            # Snakemake orchestration
│   ├── Snakefile
│   ├── rules/           # per-stage propose/execute rule pairs + manifest
│   └── envs/            # (reserved) per-rule conda env YAMLs — see envs/README.md
└── tests/               # (empty placeholder) unit tests are planned, see the approved design
```

Per-run state (artifacts, proposals, checkpoints, deliverables, internal) is written
under `run_dir` from your config — never inside the source tree.

## Environment

Implementation developed against `cell_annotation` micromamba env with pip-installed
`muon`, `scrublet`, `leidenalg`, `snakemake`, `mudata`. The plan's `workflow/envs/*.yaml`
files are the canonical production conda definitions.

**Ambient-correction R dependency (optional).** S1a calls DecontX (`celda`) or
SoupX (`SoupX`) via `Rscript`. If R / the requested package isn't installed,
S1a degrades to pass-through and records `s1a_ambient.method = "skipped_no_r"`
in `parameters.yaml`; the rest of the pipeline runs normally. To enable:

```bash
Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("celda","SoupX"))'
```

Note: SnapATAC2 function names (`pp.import_fragments`, `metrics.tsse`, `pp.add_tile_matrix`,
`pp.select_features`, `tl.spectral`, `tl.leiden`, `tl.umap`) were selected for SnapATAC2
>=2.6; verify against the installed version at execute time.
