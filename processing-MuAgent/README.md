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
  autodetect for RNA and ATAC inputs (see table below), fragments validation
  (+ tbi), and a **diagnostics-driven pairing decision**: detection (direct or
  suffix-normalized barcode overlap) is advisory; user `declare-branch paired` +
  optional `barcode_translation_path` / `cell_metadata_path` are consulted via a
  ladder before committing the workflow branch. When the ladder cannot establish
  cell-level pairing, S0 auto-downgrades `paired → separate` with the reason
  surfaced in `validation_report.json`. No barcode pre-intersection at S0 — S3
  is the sole enforcement point for the paired branch. Metadata handling
  (minimal reconstruction when absent) is unchanged.

  **Supported RNA input formats (`rna_path`):**

  | Format tag | File pattern | Notes |
  |------------|-------------|-------|
  | `10x_h5` | `*.h5` | Cell Ranger (ARC) HDF5; GEX features filtered automatically |
  | `10x_mex` | directory | 10x MEX bundle with `matrix.mtx[.gz]` + `barcodes.tsv[.gz]` |
  | `h5ad` | `*.h5ad` | AnnData; `.X` must contain raw integer counts |
  | `dense_txt` | `*.txt.gz`, `*.tsv.gz` | Dense genes × cells tab-delimited matrix (common GEO supplementary format). Row 0 = cell-barcode header; rows 1+ = gene symbol + counts. Loaded in 500-gene chunks to bound peak RAM. |

  **Supported ATAC input formats (`atac_fragments_path`):**

  | Format tag | File pattern | Notes |
  |------------|-------------|-------|
  | `fragments_tsv` | `*.tsv.gz` + `*.tsv.gz.tbi` | Standard 5-column bgzipped fragments file (`chrom start end barcode count`); tabix index must be present |
  | `bed4` *(auto-convert)* | `*.bed.gz` | 4-column BED (`chrom start end barcode`). S0 auto-converts to a standard 5-column `fragments.tsv.gz` using `zcat → awk → sort → bgzip → tabix`. The source file is **never modified**; the derived `.tsv.gz` + `.tbi` are written alongside it. Windows `\r\n` line endings in the source are handled automatically. Requires `bgzip` and `tabix` (htslib) on PATH. |
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
  (union / intersection). **Default is union** (cells flagged by either
  detector removed) for `study_goal=clustering_inference` or unspecified;
  intersection is recommended only when `study_goal=rare_populations`. The
  user confirms the policy at the S3 checkpoint. Raw per-detector calls are
  preserved in `calls.parquet`. On the paired branch, S3 also performs the
  joint barcode intersection (see "Paired multiome" below) so downstream
  stages operate on the joint cell set; a sentinel `joint_barcodes.txt` is
  written alongside the post-doublet h5ads.
- **S4 RNA norm + HVG** — log-normalize (target_sum=1e4) + HVG (`seurat_v3` on counts).
- **S5 ATAC TF-IDF + LSI and Peak Matrix Export** — Performs TF-IDF normalization and spectral embedding (LSI) on the SnapATAC2 tile matrix (`bin_size=500`, unified with S3). In parallel, the stage attempts to export a feature (cell-by-feature) matrix for ATAC using the following priority order:
  0. **User-supplied peaks:** if `atac_peaks_path` is set in `run.yaml`, S5 builds the peak-by-cell matrix from those intervals via SnapATAC2's `make_peak_matrix` ("user_peaks" mode). Highest priority — user intent wins.
  1. **ARC peak matrix:** when `single_file_multiome` was detected at S0 (Cell Ranger ARC combined `.h5`), the pre-called peaks are used directly ("arc_h5" mode).
  2. **MACS3 from fragments:** otherwise, peaks are called from fragments via SnapATAC2's MACS3 integration ("macs3_from_fragments" mode).
  3. **Tile-matrix fallback:** if all peak paths fail, the verified SnapATAC2 tile matrix is exported ("tile_matrix_fallback" mode).
  The LSI embedding (used by S6/S7/S8 for clustering and UMAP) is computed from the tile matrix regardless and is unchanged.
- **S6 Dim reduction + neighbors** — For RNA: applies `sc.pp.scale` (optional), then PCA; the number of principal components (`n_pcs`) is determined using a chord-distance ("elbow") heuristic on the explained-variance curve, capped at `rna_n_pcs_max`. Nearest-neighbors are then computed on the PCA space.
  For ATAC: cells are embedded using SnapATAC2’s LSI (from previous stage), and neighbor graphs are computed directly on the LSI representation.
- **S7 Clustering** — Leiden resolution sweep with per-modality grid, stable-region knee picker; RNA tilt=higher, ATAC tilt=lower.
- **S8 UMAP** — per-modality UMAP; paired → `processed.h5mu`, separate → two `.h5ad`.
  On the paired branch S8 expects RNA and ATAC barcodes to already match
  (S3 enforces this); the final assembly contains a defensive re-intersection
  that is logged when (and only when) it actually filters anything.
- **manifest** — `run_manifest.json` (handoff contract v1.0.0).

## Paired multiome

The paired branch admits three input shapes:

1. A single Cell Ranger ARC `.h5` (combined GEX + Peaks; barcodes share by construction).
2. Cell Ranger GEX `.h5` + ATAC fragments where the two whitelists match directly
   (or differ only by a `-N` / `_LIBRARY` suffix).
3. **Independent GEX + ATAC pipelines whose barcodes live in different 10x whitelists.**
   This requires an optional 2-column TSV at `barcode_translation_path` (or a
   `cell_metadata_path` exposing both `rna_barcode` and `atac_barcode` columns) so
   S0 can rewrite ATAC barcodes into RNA-space before any QC runs.

In all three cases, the final `processed.h5mu` contains only cells passing
both RNA and ATAC QC with matching barcodes across modalities.

**Diagnostics ladder (S0):** Detection at S0 is *advisory*; the committed
`workflow_branch` is decided by a ladder that takes the first rung that
succeeds:

1. Direct barcode overlap ≥0.99 → paired (`pairing.exact_barcode_match`).
2. Suffix-normalized overlap ≥0.99 → paired (`pairing.prefix_suffix_normalized`).
3. `barcode_translation_path` translation, then overlap ≥0.99 → paired
   (`pairing.translation_table`). Translation parquet is persisted at
   `internal/artifacts/s0_ingest/barcode_translation.parquet` and read by S2
   to produce a one-time translated copy of `atac_fragments.tsv.gz` upstream
   of the SnapATAC2 import call.
4. `cell_metadata_path` exposing both `rna_barcode` + `atac_barcode` columns
   → treated as a translation table, same rule as rung 3.
5. None of the above succeed → committed branch downgrades to `separate`
   with the reason recorded in `validation_report.json#pairing.downgrade_reason`.

If the user declared `paired` via `executor declare-branch paired` but the
ladder commits `separate`, the stage **does not crash** — the report flags the
downgrade and the agent surfaces it verbatim per `system_prompt.md` hard rule 3.

**Barcode Intersection enforcement:** No barcode pre-intersection at S0 — S1
and S2 QC each see their full modality barcode set, so QC thresholds remain
modality-native. After S3 doublet removal, the paired branch intersects RNA
and ATAC survivor sets; the joint set is written to both
`rna_post_doublet.h5ad` and `atac_post_doublet.h5ad`, and exported as
`joint_barcodes.txt`. An empty intersection at S3 raises with a remediation
message pointing at the pairing decision and the QC thresholds. Finally, at
S8 assembly, the MuData writer re-checks barcode equality *before* MuData
construction; an empty intersection is a hard error, a partial mismatch
triggers a logged subset.

**Doublet removal policy (union by default).** Doublet flagging runs independently per modality (Scrublet for RNA, SnapATAC2 scrublet for ATAC); the two flag sets are then reconciled at S3. The default is **union** — a cell is removed if flagged by *either* detector (`study_goal=clustering_inference` or unspecified). Switching `study_goal` to `rare_populations` recommends **intersection** instead (remove only cells flagged by *both* detectors), which trades a slightly higher residual-doublet rate for a much lower false-positive removal rate. The chosen policy is surfaced in `plan_review.md` and confirmed at the S3 checkpoint; per-detector scores and boolean flags are preserved in `calls.parquet`, so an alternative cut can be applied retrospectively without re-running S3.

**Diagnostic vs final clustering.** RNA-only clustering (`leiden_rna`) and
ATAC-only clustering (`leiden_atac`) at S7 are run independently per
modality and may yield different numbers of clusters — they are
**diagnostic**, not the final joint clustering. Both are computed on the
same joint cell set defined at S3, so labels can be cross-tabulated
directly. Joint clustering (e.g. WNN / MOFA+) is out of scope for this
preprocessing stage.

**Where to look in artifacts:**

- `internal/artifacts/s3_doublets/joint_barcodes.txt` — sentinel joint set
  (paired branch only).
- `internal/parameters.yaml` key `s3_doublets.paired_intersection` —
  `n_joint`, `n_dropped_rna_at_join`, `n_dropped_atac_at_join`.
- `deliverables/post_run/qc_summary.md` flow table — row "6. after S3
  doublet removal + paired intersection" shows the joint count; row 7
  shows the per-stage delta to the final object.

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
genome_assembly:       GRCh38   # or mm10

# --- RNA input (any one of the supported formats) -------------------------
rna_path:              /path/to/filtered_feature_bc_matrix.h5      # 10x h5
# rna_path:            /path/to/filtered_feature_bc_matrix/         # 10x MEX dir
# rna_path:            /path/to/counts.h5ad                         # AnnData
# rna_path:            /path/to/raw_count_matrix.txt.gz             # dense genes×cells GEO matrix

# --- ATAC input (any one of the supported formats) ------------------------
atac_fragments_path:   /path/to/atac_fragments.tsv.gz   # standard 5-col bgzip+tabix
# atac_fragments_path: /path/to/fragments.bed.gz         # 4-col BED → auto-converted by S0

# --- Optional paired-multiome inputs (all default to unset / None) --------
# 2-column TSV (rna_barcode, atac_barcode) mapping cell pairs across whitelists.
# Required only when GEX and ATAC pipelines used different 10x whitelists.
barcode_translation_path:  /path/to/barcode_translation.tsv
# Optional BED of peak intervals; consumed by S5 as the highest-priority peak
# source. If unset, S5 falls back to ARC h5 peaks (if present), then MACS3 on
# fragments, then a verified tile-matrix fallback.
atac_peaks_path:           /path/to/peaks.bed
# Optional per-cell metadata TSV (any columns; needs a `barcode` column for the
# join key). Left-joined into RNA and ATAC obs at S8 before final write. If the
# file additionally exposes `rna_barcode`+`atac_barcode` columns, S0's pairing
# ladder will treat those columns as a translation table.
cell_metadata_path:        /path/to/cell_metadata.tsv
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

## Running on HPC (PBS Pro or SLURM)

For large datasets the workflow runs in three phases. **Planning** (`p1_context`,
`p2_plan`, `plan_review`, `s0_ingest`) stays on the login node so any pairing-
detection conflict surfaces interactively before any cluster job is dispatched.
The **main preprocessing middle** (S1a–S7 propose) runs as a scheduler head-job.
Then a brief **resolution review + finish** completes S7_execute → S8 → manifest.

### One-time setup

```bash
# Imperial RDS (PBS Pro):
export PMA_PBS_QUEUE=v1_throughput72            # your allocation's queue
export PMA_PBS_PROJECT=<your project code>      # if your queue needs -P
export PMA_NOTIFY_EMAIL=<you@example.com>

# Generic SLURM site:
export PMA_SLURM_PARTITION=cpu
export PMA_SLURM_ACCOUNT=<your account>

# Optional — scale per-rule mem and walltime for larger cohorts (default 1):
export PMA_RESOURCES_SCALE=2
```

Activate the project conda env (`grn` by default; set `PMA_CONDA_ENV` to override).

### Phase A — planning (login node, inside `tmux`)

```bash
tmux new -s pma
processing-muagent init --config config/run.yaml
CFG=<run_dir>/deliverables/pre_run/config/run.yaml
processing-muagent run --config $CFG --target s0_ingest_execute
# Walk through plan_review approval / branch declaration as in the local flow.
```

### Phase B — submit the unattended head-job

```bash
processing-muagent submit --config $CFG --executor pbs \
    --auto-approve --auto-approve-except s7_clustering
```

`--auto-approve-except s7_clustering` keeps the clustering-resolution gate
honoured. The head-job stops at `s7_clustering` propose and emails
`$PMA_NOTIFY_EMAIL` (if set). Or watch from a terminal:

```bash
processing-muagent status --watch --config $CFG
```

### Phase C — review + finish

Open in any browser:

```
<run_dir>/deliverables/post_run/notebooks/resolution_review.html
```

The accompanying `resolution_review.ipynb` is for power users who want to
re-cluster at custom resolutions interactively. Approve or revise:

```bash
processing-muagent approve s7_clustering --config $CFG
# OR revise a resolution:
processing-muagent revise s7_clustering s7_clustering.rna.resolution=1.2 --config $CFG

processing-muagent submit --config $CFG --executor pbs
```

### Foreground cluster mode (alternative to `submit`)

If you'd rather keep snakemake in your tmux session (lowest-latency approvals,
no head-job queue time):

```bash
processing-muagent run --config $CFG --executor pbs
```

Snakemake stays attached to the terminal and dispatches per-rule cluster jobs;
it exits cleanly when it hits an unapproved gate. Re-invoke after each approval.

### Per-stage resources

Edit `workflow/resources.smk` if you need to override mem/walltime/cpus. The
table is the single source of truth for both PBS and SLURM profiles. OOM-killed
jobs are retried once at 2× memory (`restart-times: 1`).

## Repository layout

```
processing-MuAgent/
├── config/              # example run configurations
├── executor/            # Python implementation (stages, methods, CLI, helpers)
│   ├── stages/          # per-stage scripts S0..S8 (+ s7_notebook deliverable builder)
│   ├── methods/         # named-method helpers (MAD thresholds, resolution sweep, doublet policy)
│   └── hpc.py           # PBS/SLURM head-job submission helpers
├── workflow/            # Snakemake orchestration
│   ├── Snakefile        # localrules: declared for planning + propose + manifest
│   ├── resources.smk    # per-stage mem/runtime/cpus (single source of truth)
│   ├── rules/           # per-stage propose/execute rule pairs + manifest
│   ├── envs/            # conda env (workflow/envs/processing.yaml mirrors `grn`)
│   └── profiles/
│       ├── pbs/         # PBS Pro snakemake profile (qsub wrapper)
│       └── slurm/       # SLURM snakemake profile
├── scripts/             # launch_runner.sh + head-job templates (runner.{pbs,slurm})
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
