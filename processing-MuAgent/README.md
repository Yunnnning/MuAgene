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
  metadata handling (minimal reconstruction when absent).
- **S1 RNA QC** — MAD-derived thresholds on total_counts / n_genes / pct_counts_mt.
- **S2 ATAC QC** — TSS enrichment + nucleosome signal + fragment-count MAD via SnapATAC2.
- **S3 Doublets** — Scrublet (RNA) + ATAC doublet heuristic; four-way overlap summarised
  and goal-based removal-policy recommendation (union / intersection). Raw calls preserved.
- **S4 RNA norm + HVG** — log-normalize (target_sum=1e4) + HVG (`seurat_v3` on counts).
- **S5 ATAC TF-IDF + LSI** — SnapATAC2 tile matrix + spectral embedding, drop first component.
- **S6 Dim reduction + neighbors** — RNA PCA + neighbors; ATAC neighbors on LSI 2..N.
- **S7 Clustering** — Leiden resolution sweep `[0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2]`,
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

Note: SnapATAC2 function names (`pp.import_fragments`, `metrics.tsse`, `pp.add_tile_matrix`,
`pp.select_features`, `tl.spectral`, `tl.leiden`, `tl.umap`) were selected for SnapATAC2
>=2.6; verify against the installed version at execute time.
