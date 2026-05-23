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
  (union / intersection). **Default is union** (cells flagged by either
  detector removed) for `study_goal=clustering_inference` or unspecified;
  intersection is recommended only when `study_goal=rare_populations`. The
  user confirms the policy at the S3 checkpoint. Raw per-detector calls are
  preserved in `calls.parquet`. On the paired branch, S3 also performs the
  joint barcode intersection (see "Paired multiome" below) so downstream
  stages operate on the joint cell set; a sentinel `joint_barcodes.txt` is
  written alongside the post-doublet h5ads.
- **S4 RNA norm + HVG** — log-normalize (target_sum=1e4) + HVG (`seurat_v3` on counts).
- **S5 ATAC TF-IDF + LSI and Peak Matrix Export** — Performs TF-IDF normalization and spectral embedding (LSI) on the SnapATAC2 tile matrix (`bin_size=500`, unified with S3). In parallel, the stage attempts to export a feature (cell-by-feature) matrix for ATAC using the following approach:
  1. **Peak Matrix:** If available, loads the peak matrix directly from Cell Ranger ARC (10x) multiome `.h5` files ("arc_h5" mode). If ARC peak matrix is unavailable, generates peaks from fragments using SnapATAC2’s MACS3 integration, then constructs a peak-by-cell matrix ("macs3_from_fragments" mode).
  2. **Fallback:** If both peak matrix approaches fail, exports the verified SnapATAC2 tile matrix ("tile_matrix_fallback" mode).
  The resulting exported matrix is used for downstream integration; the LSI embedding is always preserved for ATAC clustering.
- **S6 Dim reduction + neighbors** — For RNA: applies `sc.pp.scale` (optional), then PCA; the number of principal components (`n_pcs`) is determined using a chord-distance ("elbow") heuristic on the explained-variance curve, capped at `rna_n_pcs_max`. Nearest-neighbors are then computed on the PCA space.
  For ATAC: cells are embedded using SnapATAC2’s LSI (from previous stage), and neighbor graphs are computed directly on the LSI representation.
- **S7 Clustering** — Leiden resolution sweep with per-modality grid, stable-region knee picker; RNA tilt=higher, ATAC tilt=lower.
- **S8 UMAP** — per-modality UMAP; paired → `processed.h5mu`, separate → two `.h5ad`.
  On the paired branch S8 expects RNA and ATAC barcodes to already match
  (S3 enforces this); the final assembly contains a defensive re-intersection
  that is logged when (and only when) it actually filters anything.
- **manifest** — `run_manifest.json` (handoff contract v1.0.0).

## Paired multiome

The paired branch (single Cell Ranger ARC `.h5` or matched RNA matrix +
ATAC fragments) guarantees that **the final object contains only cells
passing both RNA and ATAC QC, with matching barcodes across modalities**.

**Barcode Intersection:** In the paired workflow branch, barcode alignment is enforced in three steps to ensure that only cells passing both RNA and ATAC QC are included. First, during S0 ingest, RNA cells are filtered to those with barcodes present in the ATAC fragments, while the ATAC set is left unfiltered to maximize statistical power for ATAC QC. Next, after S1 and S2 QC and S3 doublet removal, a strict intersection of surviving RNA and ATAC barcodes is performed; the resulting joint set is written to both `rna_post_doublet.h5ad` and `atac_post_doublet.h5ad`, and exported as `joint_barcodes.txt` for downstream validation. Finally, at S8 assembly, the MuData writer re-checks that the RNA and ATAC barcodes are identical; if any mismatch arises due to cell losses in earlier stages, the intersection is recalculated and the event is logged to maintain consistency.

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
