# processing-MuAgent

Multiome (scRNA-seq + scATAC-seq) preprocessing subagent. Takes raw 10x Genomics
multiome outputs and produces QC'd, dim-reduced, clustered, UMAP-visualised data
**per modality**, then **stops**.

## Hard stop (out of scope)

No integration, no WNN, no label transfer, no cell annotation, no marker genes,
no GRN. Those belong to downstream subagents.

## Workflow

Pre-preprocessing phases (run before any QC):

- **P1 Context Extraction** — Biological Context Report (organism, tissue, assay, DOIs)
  + DOI-based prior-analysis extraction. Conflicts surfaced, never auto-resolved.
- **P2 Preprocessing Plan Generation** — holistic `preprocessing_plan.json` for every
  downstream stage; approved atomically before execution.

Preprocessing stages:

- **S0 Ingest** — format autodetect (10x h5 / MEX / h5ad / custom), fragments
  validation (+ tbi), pairing detection (paired vs separate workflow branches),
  metadata handling (minimal reconstruction when absent). Accepts both
  Cell Ranger **filtered** and **raw** matrices; raw matrices trigger a
  barcode-rank knee for cell calling. Optional `rna_raw_path` carries the
  raw drops alongside the filtered cells (used by SoupX in S1a).
- **S1a Ambient RNA correction** — DecontX (filtered counts only) or SoupX
  (raw + filtered) auto-dispatched from `s0` outputs. Pass-through when
  R / Bioconductor is unavailable; deviation explicitly recorded in provenance.
  Per-cell rho is exposed in `obs['ambient_contamination']` and
  `s1a_ambient/contamination.parquet`; raw counts are preserved in
  `layers['counts_raw']`.
- **S1 RNA QC** — MAD-derived thresholds on total_counts / n_genes / pct_counts_mt
  + `pct_counts_ribo` ceiling, computed on the decontaminated counts from S1a.
- **S2 ATAC QC** — TSS enrichment + per-cell nucleosome signal (Signac-style
  `mono / nucleosome_free`) + fragment-count MAD via SnapATAC2.
- **S3 Doublets** — Scrublet (RNA, sparse-CSR input, adaptive
  `expected_doublet_rate ≈ 0.0008 × n_cells`) + SnapATAC2 scrublet (ATAC);
  four-way overlap summarised and goal-based removal-policy recommendation
  (union / intersection). Raw calls preserved.
- **S4 RNA norm + HVG** — log-normalize (target_sum=1e4) + HVG (`seurat_v3` on counts).
- **S5 ATAC TF-IDF + LSI** — SnapATAC2 tile matrix (`bin_size=500`, unified
  with S3) + spectral embedding, drop first component.
- **S6 Dim reduction + neighbors** — RNA `sc.pp.scale` + PCA with chord-distance
  elbow-detected `n_pcs` (capped at `rna_n_pcs_max`) + neighbors; ATAC neighbors on LSI 2..N.
- **S7 Clustering** — Leiden resolution sweep with per-modality grid,
  stable-region knee picker; RNA tilt=higher, ATAC tilt=lower.
- **S8 UMAP** — per-modality UMAP; paired → `processed.h5mu`, separate → two `.h5ad`.
- **manifest** — `run_manifest.json` (handoff contract v1.0.0).

## Execution engine

Snakemake: each stage is a `<stage>_propose` + `<stage>_execute` rule pair, with
`<stage>_execute` gated on `checkpoints/<stage>.approved`.

## CLI

```bash
# one-time: install the executor as an importable package
pip install -e .

# scaffold a run
processing-muagent init --config config/run.yaml

# edit <run_dir>/biological_context.md (optional) and run it
processing-muagent run --config <run_dir>/run.yaml --auto-approve
```

Interactive mode:

```bash
processing-muagent propose p1_context --config <run_dir>/run.yaml
# review proposals/p1_context.yaml
processing-muagent approve p1_context --config <run_dir>/run.yaml
processing-muagent propose s0_ingest --config <run_dir>/run.yaml
# ... etc.
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
