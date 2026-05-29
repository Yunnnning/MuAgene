# processing-MuAgent

Multiome (scRNA-seq + scATAC-seq) preprocessing subagent. Takes filtered or raw 10x Genomics multiome outputs and performs QC, dimensionality reduction, clustering, and UMAP **per modality**, then **stops** before integration.

Supported workflow branches: `paired`, `separate`, `rna_only`, `atac_only`. Declare the branch up front with `processing-muagent declare-branch`.

## Pipeline overview

**Stage order** (Snakemake DAG — ingest validation must finish before preprocessing plan assembly, because P2 reads `internal/artifacts/s0_ingest/validation_report.json`):

```
P1 context extraction → S0 ingest validation → P2 preprocessing plan → (CHECKPOINT 1) plan_review
  → S1a ambient RNA correction → S1 RNA QC → S2 ATAC QC → S3 doublets → (CHECKPOINT 2) post_qc_review
  → S4 RNA normalization + HVG → S5 ATAC TF-IDF + LSI → S6 dimensionality reduction + neighbors 
  → S7 clustering + (CHECKPOINT 3) resolution_review → S8 UMAP → manifest
```

Each stage is a Snakemake `<stage>_propose` + `<stage>_execute` pair (except `post_qc_review`, which is propose-only). Execute rules run only after `internal/checkpoints/<stage>.approved` is written by `processing-muagent approve <stage>`.

### User checkpoints (3)

Three deliberate pauses where you review deliverables and decide before heavy downstream work continues. All other stages (`s0_ingest`, `p2_plan`, `s1a`–`S3`, `S4`–`S6`, `S8`) are normally auto-approved on HPC.

| # | Checkpoint | Gate stage | When | What you decide |
|---|------------|------------|------|-----------------|
| **1** | **Plan review** | `plan_review` | After S0 + P2, before S1 | Approve the preprocessing plan (`pre_run/summary/plan_review.md`) |
| **2** | **QC review** | `post_qc_review` | After S3, before S4/S5 | Inspect QC figures + `checkpoint/qc_review/qc_summary.md`; revise **S1/S2 thresholds** and re-run if needed; on **paired** multiome, also confirm the **S3 cross-modal doublet removal policy** (union vs intersection) documented in the summary |
| **3** | **Clustering resolution review** | `s7_clustering` | After S6, before S8 | Choose Leiden resolution per modality from sweep metrics (`checkpoint/resolution_review/`). **Separate / single-modality:** sets **final** cluster labels in processed outputs. **Paired:** **diagnostic** per-modality labels for UMAP only (not joint embedding) |

**Snakemake approval gates:** 14 stages require an internal `*.approved` sentinel (including the three checkpoints above). On HPC, all gates except checkpoints **#2** and **#3** are normally auto-approved so the batch job runs unattended between your reviews.

## Workflow stages

### Planning (pre-QC)

- **P1 Context extraction** — Biological Context Report (organism, tissue, assay, DOIs) plus DOI-based prior-analysis extraction.
- **S0 Ingest** — Accepts Cell Ranger **filtered** and **raw** matrices, auto-detecting RNA and ATAC formats (see tables below) and validating fragments files. Performs a **diagnostic barcode check for paired multiome**: S0 checks for direct barcode matches, then tries matching after removing suffixes. If those don't match, it looks for a `barcode_translation_path` or `cell_metadata_path` provided by the user. No barcode intersection is performed here. If pairing can’t be confirmed, S0 downgrades the workflow from `paired` to `separate` and logs the reason in `validation_report.json`. S0 runs after P1 and before P2; its validation report feeds into the preprocessing plan. By default, S0 runs on the login node, but will auto-retry as a **cluster job** for large datasets using the same HPC setup as later stages.

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
  | `bed4` *(auto-convert)* | `*.bed.gz` | 4-column BED (`chrom start end barcode`). S0 auto-converts to a standard 5-column `fragments.tsv.gz` using `zcat → awk → sort → bgzip → tabix`. The source file is **never modified**; the derived `.tsv.gz` + `.tbi` are written alongside it. Windows `\r\n` line endings are handled automatically. Requires `bgzip` and `tabix` (htslib) on PATH. |

- **P2 Preprocessing plan generation** — Creates `preprocessing_plan.json` with execution and parameter settings for all downstream stages, using outputs from P1 context and S0 ingest.
- **plan_review** — Generates a summary at `deliverables/pre_run/summary/plan_review.md` for the user to review. The workflow pauses here until approval, before any S1–S8 execute rule runs.

### Preprocessing

- **S1a Ambient RNA correction** — DecontX (filtered counts only) or SoupX (raw + filtered) auto-dispatched from S0 outputs. Pass-through when R / Bioconductor is unavailable.
- **S1 RNA QC** — MAD-derived thresholds on `total_counts` / `n_genes` / `pct_counts_mt` plus a `pct_counts_ribo` ceiling, computed on decontaminated counts from S1a. Writes pre/post QC violin figures to `deliverables/checkpoint/qc_review/`.
- **S2 ATAC QC** — TSS enrichment, per-cell nucleosome signal (Signac-style `mono / nucleosome_free`), and fragment-count MAD via SnapATAC2. Writes fragment-size distribution figures to `deliverables/checkpoint/qc_review/`.
- **S3 Doublets** — Per-modality doublet detection, then branch-specific reconciliation:
  - **RNA:** Scrublet (sparse-CSR input; `expected_doublet_rate ≈ 0.0008 × n_cells`, capped at 10%).
  - **ATAC:** SnapATAC2 scrublet (thresholds configurable in the preprocessing plan).
  - **separate / single-modality branches:** Each modality is filtered independently by its own detector; per-modality calls are saved in `calls.parquet`.
  - **paired branch:** Also performs joint barcode intersection after doublet removal; the applied cross-modal policy (union vs intersection) is confirmed at the **QC review checkpoint** (`checkpoint/qc_review/qc_summary.md`).
- **post_qc_review** — **QC review checkpoint (#2).** Propose-only gate between S3 and dimensionality reduction. Generates doublet histograms, cell-count waterfall, and `checkpoint/qc_review/qc_summary.md` (S1–S3 metrics + paired S3 policy). Revise S1/S2 thresholds or `s3_doublets.removal_policy` and re-run affected stages before approving.
- **S4 RNA norm + HVG** — Log-normalize (`target_sum=1e4`) + HVG selection (`seurat_v3` on counts).
- **S5 ATAC TF-IDF + LSI and peak matrix export** — TF-IDF normalization and spectral embedding (LSI) on the SnapATAC2 tile matrix (`bin_size=500`, unified with S3). In parallel, exports a feature (cell-by-feature) matrix using this priority order:
  0. **User-supplied peaks** — `atac_peaks_path` in `run.yaml` → SnapATAC2 `make_peak_matrix` (`user_peaks` mode).
  1. **ARC peak matrix** — pre-called peaks from a combined Cell Ranger ARC `.h5` detected at S0 (`arc_h5` mode).
  2. **MACS3 from fragments** — SnapATAC2 MACS3 integration (`macs3_from_fragments` mode).
  3. **Tile-matrix fallback** — verified SnapATAC2 tile matrix (`tile_matrix_fallback` mode).

  LSI embedding (used by S6–S8) is always computed from the tile matrix regardless of peak-export mode.
- **S6 Dim reduction + neighbors** — **RNA:** optional `sc.pp.scale`, then PCA; `n_pcs` from a chord-distance elbow on explained variance, capped at `rna_n_pcs_max`; nearest-neighbors on PCA space. **ATAC:** SnapATAC2 LSI embedding from S5; neighbor graph on LSI.
- **S7 Clustering** — Leiden resolution sweep with per-modality grid and stable-region knee picker. **Resolution review checkpoint (#3):** `checkpoint/resolution_review/resolution_review.html` / `.ipynb`. Separate branch: chosen resolutions become final labels. Paired branch: diagnostic per-modality labels for UMAP only.
- **S8 UMAP** — Per-modality UMAP. **Paired** → `processed.h5mu`; **separate** → `rna_processed.h5ad` + `atac_processed.h5ad`. On the paired branch, S8 expects matching barcodes from S3; final assembly includes a defensive re-intersection logged only when it filters cells.
- **manifest** — `run_manifest.json` handoff contract (v1.0.0), final `qc_summary.md`, and `layout.json`.

## Paired multiome

The paired branch admits three input shapes:

1. A single Cell Ranger ARC `.h5` (combined GEX + Peaks; barcodes match by construction).
2. Cell Ranger GEX `.h5` + ATAC fragments where cell barcodes match directly (or differ only by a `-N` / `_LIBRARY` suffix).
3. Independent GEX + ATAC pipelines whose barcodes live in different 10x whitelists — requires a 2-column TSV at `barcode_translation_path` (or `cell_metadata_path` with `rna_barcode` + `atac_barcode` columns) so S0 can rewrite ATAC barcodes into RNA space before QC.

In all cases, the final `processed.h5mu` contains only cells passing both RNA and ATAC QC with matching barcodes.

### Barcode matching and intersection

In paired RNA and ATAC workflows, the pipeline first checks how well the barcodes from each dataset align:

1. **Strong match:** If at least 80% of barcodes directly overlap between RNA and ATAC, the data is treated as paired.
2. **Subset match:** If most ATAC barcodes are found in the RNA set (or vice versa), the data is also paired.
3. **Minor differences:** If barcodes are similar but differ by prefixes or suffixes, they are normalized and re-checked. If 80% still match, pairing proceeds.
4. **Translation table:** If barcodes are different but you provide a translation table (`barcode_translation_path` or `cell_metadata_path`), this is used to match barcodes.
5. **Ambiguous match:** If the overlap is between 30% and 80% and can't be resolved above, you must declare if the data is paired or provide a translation.
6. **Low match:** If none of the above, the data is treated as separate.

If you select paired mode but the automatic check does not support it, you will be notified to review the findings; the process will not stop, but a report will flag the issue.

**Barcode intersection enforcement:**

- **S0:** Barcodes are checked but not intersected. S1 and S2 each see the full modality barcode set in RNA and ATAC modalities.
- **S3 (paired):** After doublet removal, only barcodes found in both modalities are kept. Empty intersection raises with a remediation message.
- **S8 (paired):** Before final output, the workflow re-checks barcode equality before construction; empty intersection is a hard error; partial mismatch triggers a logged subset.

**Doublet removal in paired data:**

Doublets are detected separately in RNA and ATAC; their results are combined during S3 reconciliation. By default (**union mode**), a cell is removed if either detector (RNA or ATAC) flags it as a doublet (`study_goal=clustering_inference` or unset). For `study_goal=rare_populations`, **intersection mode** is recommended: remove only if both detectors flag the cell. Detector scores and flags are saved in `calls.parquet` for later review.

### Diagnostic vs final clustering

At S7, RNA-only (`leiden_rna`) and ATAC-only (`leiden_atac`) clustering are run separately for diagnostic comparison—not for joint clustering. In paired mode, both clusterings use the same intersected cell set from S3; in separate mode, they use their respective modality’s cells. Joint multimodal clustering (e.g., WNN or MOFA+) is not performed in this preprocessing workflow.

## Repository layout

```
processing-MuAgent/
├── agent/               # chat-runtime prompts (system_prompt, interaction_flow)
├── config/              # example run configurations
├── executor/            # Python implementation (stages, methods, CLI, helpers)
│   ├── stages/          # per-stage scripts S0..S8 + post_qc_review
│   ├── methods/         # MAD thresholds, resolution sweep, doublet policy
│   └── hpc.py           # PBS/SLURM head-job submission helpers
├── workflow/            # Snakemake orchestration
│   ├── Snakefile        # localrules for planning + propose + manifest
│   ├── resources.smk    # per-stage mem/runtime/cpus
│   ├── rules/           # per-stage propose/execute rule pairs + manifest
│   ├── envs/            # conda env (mirrors `grn`)
│   └── profiles/
│       ├── pbs/         # PBS Pro snakemake profile
│       └── slurm/       # SLURM snakemake profile
└── scripts/             # launch_runner.sh + head-job templates
```


## Running on HPC (PBS Pro or SLURM)

On a cluster, heavy compute stages run as scheduler jobs (PBS Pro or SLURM). The agent drives the workflow through the same checkpoints as local mode; you only need to configure your site once (below). Everything else — init, submit, approve, revise — is handled via the CLI or the chat agent.

### One-time setup

Set these environment variables for your cluster before the first run:

```bash
# PBS Pro example:
export PMA_PBS_QUEUE=<your_queue_name>
export PMA_PBS_PROJECT=<your_project_code>
export PMA_NOTIFY_EMAIL=<your_email_address>

# SLURM example:
export PMA_SLURM_PARTITION=<your_partition_name>
export PMA_SLURM_ACCOUNT=<your_account_name>

# Optional — scale per-rule memory and walltime (default is 1):
export PMA_RESOURCES_SCALE=2
```

`PMA_NOTIFY_EMAIL` is optional but recommended: you receive mail when a submitted batch finishes or pauses at a review checkpoint. For larger datasets, increase `PMA_RESOURCES_SCALE` (e.g. `2` for ~30k cells, `4` for ~100k). Per-stage CPU, memory, and walltime defaults live in `workflow/resources.smk`; OOM-killed jobs are retried once at double memory.

### Requirements

- **SLURM:** Requires Snakemake version 9 or higher. PBS Pro does not have this requirement.

### How the HPC run proceeds

| Step | Stages | Executes on | You |
|------|--------|-------------|-----|
| Planning | P1 → P2 | Login node (default), or `srun` on a compute node if the login node memory is limited| — |
| S0 ingest | S0 | Login node (default); Cluster if local fails for large dataset | — |
| Checkpoint **#1** | plan_review | Login node | Review plan |
| QC | S1a → S1 → S2 → S3 | Cluster | — |
| Checkpoint **#2** | post_qc_review | — | Review QC |
| Dimred + clustering | S4 → S5 → S6 → S7 (sweep) | Cluster | — |
| Checkpoint **#3** | s7_clustering | — | Review resolution |
| Finish | S7 (labels) → S8 → manifest | Cluster | — |

**Flexible S0:** the agent runs ingest locally first. On OOM or walltime failure it configures HPC (if needed), sources `hpc.env`, and retries `s0_ingest_execute` on the cluster before continuing with P2.

Each heavy `_execute` stage runs as its own scheduler job. Gates between your reviews are auto-approved; email notification fires when a batch pauses or completes (if `PMA_NOTIFY_EMAIL` is set).

### Submit workflow

After checkpoint **#1** (`plan_review`), source `deliverables/pre_run/config/hpc.env`, then use `processing-muagent submit` (not `run`) to dispatch the Snakemake head-job. **`submit` auto-infers the Snakemake target** from checkpoint state — you do not need to pick `post_qc_review_propose`, `s7_clustering_propose`, or `all` manually. After each approval, run `submit` again and it stops at the next gate:

| Checkpoint state | Inferred target | Runs through |
|------------------|-----------------|--------------|
| `post_qc_review` not approved | `post_qc_review_propose` | S1a → S3 + QC summary, then pauses |
| `s7_clustering` not approved | `s7_clustering_propose` | S4 → S6 + resolution sweep, then pauses |
| Both approved | `all` | S7 labels → S8 → manifest |

Override with `--target <name>` only when debugging.

```bash
CFG=<run_dir>/deliverables/pre_run/config/run.yaml
source <run_dir>/deliverables/pre_run/config/hpc.env

# First heavy batch after plan review (honours QC + resolution checkpoints):
processing-muagent submit --config $CFG --executor slurm \
  --auto-approve --auto-approve-except post_qc_review \
  --auto-approve-except s7_clustering

# After QC review:
processing-muagent approve post_qc_review --config $CFG
processing-muagent submit --config $CFG --executor slurm \
  --auto-approve --auto-approve-except s7_clustering

# After resolution review:
processing-muagent approve s7_clustering --config $CFG
processing-muagent submit --config $CFG --executor slurm
```

Poll with `processing-muagent status --watch --config $CFG`. `--auto-approve-except` syntax is unchanged; repeat the flag for each gate you want to keep interactive.


## Run directory layout

Per-run state lives under `run_dir` from your config — never inside the source tree.

```
<run_dir>/
  deliverables/
    pre_run/
      config/
        run.yaml                  ← canonical config (use this for all CLI calls)
        biological_context.md     ← Biological Context Report
      summary/
        context_summary.md        ← P1 output
        plan_summary.md           ← P2 output
        plan_review.md            ← plan review gate
    checkpoint/
      qc_review/                  ← QC review checkpoint (#2): figures + qc_summary.md
      resolution_review/          ← resolution_summary.md + resolution_review.{html,ipynb}
    post_run/                     ← flat final deliverables
      s8_umap_*.{png,pdf}         ← UMAP figures only
      processed.h5mu              ← or rna/atac_processed.h5ad (separate branch)
      review_processed_h5mu.{ipynb,py}
      qc_summary.md               ← final QC summary
      run_manifest.json           ← handoff artifact
      layout.json
  internal/
    artifacts/sN_<stage>/         ← intermediate stage outputs
    proposals/                    ← <stage>.yaml + awaiting_approval sentinels
    checkpoints/                  ← <stage>.approved sentinels
    parameters.yaml
    state.yaml
    log.jsonl
```

## CLI

**Step 1 — Install** (from inside `processing-MuAgent/`):

```bash
cd /path/to/processing-MuAgent
pip install -e .
```

**Step 2 — Edit the example config** at `config/run.example.yaml`. At minimum:

```yaml
run_dir:               /path/to/your/output/run_01
genome_assembly:       GRCh38   # or mm10
study_goal:            clustering_inference   # or rare_populations

# --- RNA input (any supported format) -----------------------------------
rna_path:              /path/to/filtered_feature_bc_matrix.h5
# rna_raw_path:        /path/to/raw_feature_bc_matrix.h5   # enables SoupX in S1a

# --- ATAC input ---------------------------------------------------------
atac_fragments_path:   /path/to/atac_fragments.tsv.gz
# atac_fragments_path: /path/to/fragments.bed.gz            # auto-converted by S0

# --- Optional paired-multiome inputs ------------------------------------
barcode_translation_path:  /path/to/barcode_translation.tsv   # rna_barcode, atac_barcode
atac_peaks_path:           /path/to/peaks.bed                 # highest-priority peak source for S5
cell_metadata_path:        /path/to/cell_metadata.tsv         # obs join at S8; pairing ladder if rna+atac cols
```

**Step 3 — Scaffold the run directory:**

```bash
processing-muagent init --config config/run.example.yaml
```

`init` creates `<run_dir>/` and copies your config to `<run_dir>/deliverables/pre_run/config/run.yaml`. It also writes the Biological Context Report template at `deliverables/pre_run/config/biological_context.md`.

**Step 4 — Declare branch and run:**

```bash
CFG=<run_dir>/deliverables/pre_run/config/run.yaml

processing-muagent declare-branch paired --config $CFG   # or separate | rna_only | atac_only

# Fully automatic (honours all gates unless you exclude them):
processing-muagent run --config $CFG --auto-approve

# Check status at any point:
processing-muagent status --config $CFG
```

**Interactive / checkpoint-by-checkpoint mode:**

```bash
CONFIG=<run_dir>/deliverables/pre_run/config/run.yaml

processing-muagent propose p1_context --config $CONFIG
# review: <run_dir>/internal/proposals/p1_context.yaml
processing-muagent approve p1_context --config $CONFIG

processing-muagent propose s0_ingest --config $CONFIG
# review: <run_dir>/internal/proposals/s0_ingest.yaml
#         <run_dir>/internal/artifacts/s0_ingest/validation_report.json
processing-muagent approve s0_ingest --config $CONFIG

processing-muagent propose p2_plan --config $CONFIG
# review: <run_dir>/deliverables/pre_run/summary/plan_summary.md
processing-muagent approve p2_plan --config $CONFIG

processing-muagent plan-review --config $CONFIG
# review: <run_dir>/deliverables/pre_run/summary/plan_review.md
processing-muagent approve plan_review --config $CONFIG

# S1a → S8:
for STAGE in s1a_ambient s1_rna_qc s2_atac_qc \
             s3_doublets post_qc_review s4_rna_norm s5_atac_lsi \
             s6_dimred s7_clustering s8_umap; do
    processing-muagent propose $STAGE --config $CONFIG
    # review: <run_dir>/internal/proposals/$STAGE.yaml
    # post_qc_review (QC review #2): deliverables/checkpoint/qc_review/qc_summary.md
    # s7_clustering (resolution review #3): deliverables/checkpoint/resolution_review/resolution_review.html
    processing-muagent approve $STAGE --config $CONFIG
done
```

Other useful commands:

```bash
processing-muagent revise s7_clustering s7_clustering.rna.resolution=1.2 --config $CFG
processing-muagent resolution-compare --config $CFG --rna 1.0,1.2 --atac 0.6,0.8
processing-muagent run --config $CFG --no-context   # explicit opt-out of biological context
processing-muagent hpc-info                         # probe queues/partitions on login node
```

On HPC after plan review, use `submit` (see **Submit workflow** above) instead of `run`. `submit` auto-infers the phase target; `--auto-approve-except post_qc_review --auto-approve-except s7_clustering` keeps the two QC/resolution checkpoints interactive.


## Environment

Recreate the project conda env from `workflow/envs/processing.yaml` (default name `grn`; override with `PMA_CONDA_ENV`):

```bash
micromamba env create -n grn -f workflow/envs/processing.yaml
micromamba activate grn
pip install -e .
```

**Ambient-correction R dependency (optional).** S1a calls DecontX (`celda`) or SoupX (`SoupX`) via `Rscript`. If R / the requested package is not installed, S1a degrades to pass-through and records `s1a_ambient.method = "skipped_no_r"` in `parameters.yaml`. To enable:

```bash
Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("celda","SoupX"))'
```

SnapATAC2 function names (`pp.import_fragments`, `metrics.tsse`, `pp.add_tile_matrix`, `pp.select_features`, `tl.spectral`, `tl.leiden`, `tl.umap`) were selected for SnapATAC2 ≥ 2.6; verify against the installed version at execute time.
