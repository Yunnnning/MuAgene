# processing-MuAgent

Multiome (scRNA-seq + scATAC-seq) preprocessing subagent. Takes filtered or raw 10x Genomics multiome outputs and performs QC, dimensionality reduction, clustering, and UMAP **per modality**, then **stops** before integration.

Supported workflow branches: `paired`, `separate`, `rna_only`, `atac_only`. Declare the branch up front with `processing-muagent declare-branch`.

## Pipeline overview

**Stage order** (Snakemake DAG — ingest validation must finish before preprocessing plan assembly, because P2 reads `internal/artifacts/s0_ingest/validation_report.json`):

```
P1 context extraction → S0 ingest validation → P2 preprocessing plan → plan_review
  → S1a ambient RNA correction → S1 RNA QC → S2 ATAC QC → S3 doublets → post_qc_review
  → S4 RNA normalization + HVG → S5 ATAC TF-IDF + LSI
  → S6 dimensionality reduction + neighbors → S7 clustering → S8 UMAP → manifest
              ↑                              ↑                        ↑
         CHECKPOINT 1                   CHECKPOINT 2            CHECKPOINT 3
        (plan review)                  (QC review)      (resolution review)
```

| Stage ID | Full name |
|----------|-----------|
| `p1_context` | P1 context extraction |
| `s0_ingest` | S0 ingest validation |
| `p2_plan` | P2 preprocessing plan |
| `plan_review` | Plan review (checkpoint 1) |
| `s1a_ambient` | S1a ambient RNA correction |
| `s1_rna_qc` | S1 RNA QC |
| `s2_atac_qc` | S2 ATAC QC |
| `s3_doublets` | S3 doublets |
| `post_qc_review` | Post-QC review (checkpoint 2) |
| `s4_rna_norm` | S4 RNA normalization + HVG |
| `s5_atac_lsi` | S5 ATAC TF-IDF + LSI |
| `s6_dimred` | S6 dimensionality reduction + neighbors |
| `s7_clustering` | S7 clustering (checkpoint 3) |
| `s8_umap` | S8 UMAP |
| `manifest` | Run manifest + final summaries |

Each stage is a Snakemake `<stage>_propose` + `<stage>_execute` pair (except `post_qc_review`, which is propose-only). Execute rules run only after `internal/checkpoints/<stage>.approved` is written by `processing-muagent approve <stage>`.

### User checkpoints (3)

Three deliberate pauses where you review deliverables and decide before heavy downstream work continues. All other stages (`s0_ingest`, `p2_plan`, `s1a`–`S3`, `S4`–`S6`, `S8`) are normally auto-approved on HPC.

| # | Checkpoint | Gate stage | When | What you decide |
|---|------------|------------|------|-----------------|
| **1** | **Plan review** | `plan_review` | After S0 + P2, before S1 | Approve the preprocessing plan (`pre_run/summary/plan_review.md`) |
| **2** | **QC review** | `post_qc_review` | After S3, before S4/S5 | Inspect QC figures + `checkpoint/qc_review/qc_summary.md`; revise S1/S2 thresholds and re-run if needed; on **paired** multiome, also confirm the **S3 cross-modal doublet policy** (union vs intersection) documented in the summary |
| **3** | **Clustering resolution review** | `s7_clustering` | After S6, before S8 | Choose Leiden resolution per modality from sweep metrics (`checkpoint/resolution_review/`). **Separate / single-modality:** sets **final** cluster labels in processed outputs. **Paired:** **diagnostic** per-modality labels for UMAP only (not joint embedding) |

**QC review and S3 policy:** S3 runs before the QC review checkpoint. On `paired`, the applied doublet policy and joint barcode counts are written into `qc_summary.md` — there is no separate S3 user gate. On `separate` / `rna_only` / `atac_only`, doublets are removed independently per modality; no cross-modal policy applies.

**Snakemake approval gates:** 14 stages require an internal `*.approved` sentinel (including the three checkpoints above). With `--auto-approve`, all 14 are pre-seeded; typical HPC usage keeps only checkpoints **2** and **3** open:

```bash
processing-muagent submit --config $CFG --executor pbs \
    --auto-approve --auto-approve-except post_qc_review \
    --auto-approve-except s7_clustering
```

## Workflow stages

### Planning (pre-QC)

- **P1 Context extraction** — Biological Context Report (organism, tissue, assay, DOIs) plus DOI-based prior-analysis extraction.
- **S0 Ingest** — Accepts Cell Ranger **filtered** and **raw** matrices. Format autodetect for RNA and ATAC inputs (see tables below), fragments validation (+ `.tbi`), and a **diagnostics-driven pairing decision**: detection (direct or suffix-normalized barcode overlap) is advisory; user `declare-branch` plus optional `barcode_translation_path` / `cell_metadata_path` are consulted via a ladder before committing the workflow branch. When the ladder cannot establish cell-level pairing, S0 auto-downgrades `paired → separate` with the reason in `validation_report.json`. No barcode pre-intersection at S0 — S3 is the sole enforcement point for the paired branch. **Runs after P1 and before P2** — the preprocessing plan is built from S0's validation report.

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

- **P2 Preprocessing plan generation** — Holistic `preprocessing_plan.json` for every downstream stage; assembled from P1 context + S0 ingest outputs.
- **plan_review** — Merged plan review at `deliverables/pre_run/summary/plan_review.md`: concise decision summary (8 items) plus a full parameter appendix. Hard gate before any S1–S8 execute rule runs.

### Preprocessing

- **S1a Ambient RNA correction** — DecontX (filtered counts only) or SoupX (raw + filtered) auto-dispatched from S0 outputs. Pass-through when R / Bioconductor is unavailable.
- **S1 RNA QC** — MAD-derived thresholds on `total_counts` / `n_genes` / `pct_counts_mt` plus a `pct_counts_ribo` ceiling, computed on decontaminated counts from S1a. Writes pre/post QC violin figures to `deliverables/checkpoint/qc_review/`.
- **S2 ATAC QC** — TSS enrichment, per-cell nucleosome signal (Signac-style `mono / nucleosome_free`), and fragment-count MAD via SnapATAC2. Writes fragment-size distribution figures to `deliverables/checkpoint/qc_review/`.
- **S3 Doublets** — Per-modality doublet detection, then branch-specific reconciliation:
  - **RNA:** Scrublet (sparse-CSR input; `expected_doublet_rate ≈ 0.0008 × n_cells`, capped at 10%).
  - **ATAC:** SnapATAC2 scrublet (thresholds configurable in the preprocessing plan).
  - **`separate` / single-modality branches:** Each modality is filtered independently by its own detector; per-modality calls are saved in `calls.parquet`.
  - **`paired` branch:** Also performs joint barcode intersection after doublet removal; the applied cross-modal policy (union vs intersection) is confirmed at the **QC review checkpoint** (`checkpoint/qc_review/qc_summary.md`), not at a separate S3 gate.
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
2. Cell Ranger GEX `.h5` + ATAC fragments where whitelists match directly (or differ only by a `-N` / `_LIBRARY` suffix).
3. **Independent GEX + ATAC pipelines** whose barcodes live in different 10x whitelists — requires a 2-column TSV at `barcode_translation_path` (or `cell_metadata_path` with `rna_barcode` + `atac_barcode` columns) so S0 can rewrite ATAC barcodes into RNA space before QC.

In all cases, the final `processed.h5mu` contains only cells passing both RNA and ATAC QC with matching barcodes.

### Diagnostics ladder (S0)

Detection at S0 is advisory; the committed `workflow_branch` follows the first successful rung:

1. Direct barcode overlap ≥ 0.99 → paired (`pairing.exact_barcode_match`).
2. Suffix-normalized overlap ≥ 0.99 → paired (`pairing.prefix_suffix_normalized`).
3. `barcode_translation_path` translation, then overlap ≥ 0.99 → paired (`pairing.translation_table`). Translation parquet is persisted at `internal/artifacts/s0_ingest/barcode_translation.parquet`; S2 reads it to produce a one-time translated copy of `atac_fragments.tsv.gz` before SnapATAC2 import.
4. `cell_metadata_path` with `rna_barcode` + `atac_barcode` columns → same rule as rung 3.
5. None succeed → branch downgrades to `separate`; reason in `validation_report.json#pairing.downgrade_reason`.

If the user declared `paired` but the ladder commits `separate`, the stage does not crash — the report flags the downgrade for review.

### Barcode intersection enforcement

- **S0:** No pre-intersection; S1 and S2 each see their full modality barcode set.
- **S3 (paired):** After doublet removal, intersect RNA and ATAC survivor sets. Joint set is written to both `rna_post_doublet.h5ad` and `atac_post_doublet.h5ad`, plus `joint_barcodes.txt`. Empty intersection raises with a remediation message.
- **S8 (paired):** MuData writer re-checks barcode equality before construction; empty intersection is a hard error; partial mismatch triggers a logged subset.

### Doublet removal policy (paired branch)

Doublet flagging runs independently per modality; flag sets are reconciled at S3. Default is **union** — remove if flagged by either detector (`study_goal=clustering_inference` or unspecified). With `study_goal=rare_populations`, the recommendation switches to **intersection** (remove only if flagged by both). Per-detector scores and flags are preserved in `calls.parquet` for retrospective re-cutting.

### Diagnostic vs final clustering

RNA-only (`leiden_rna`) and ATAC-only (`leiden_atac`) clustering at S7 are **diagnostic**, not joint clustering. Both run on the same joint cell set from S3 (paired) or per-modality sets (separate). Joint clustering (WNN / MOFA+) is out of scope.

### Key artifact paths

| Path | Contents |
|------|----------|
| `internal/artifacts/s3_doublets/joint_barcodes.txt` | Sentinel joint set (paired only) |
| `internal/parameters.yaml` → `s3_doublets.paired_intersection` | `n_joint`, `n_dropped_rna_at_join`, `n_dropped_atac_at_join` |
| `deliverables/checkpoint/qc_review/qc_summary.md` | QC review checkpoint — S1–S3 metrics + paired S3 doublet policy |
| `deliverables/post_run/qc_summary.md` | Final QC summary (written at manifest) |

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
        plan_review.md            ← plan review gate (summary + appendix)
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
processing-muagent approve p2_plan --config $CONFIG

processing-muagent plan-review --config $CONFIG
# review: <run_dir>/deliverables/pre_run/summary/plan_review.md  (summary + appendix)
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
```

## Running on HPC (PBS Pro or SLURM)

For large datasets the workflow runs in four phases. Planning (`p1_context`, `s0_ingest`, `p2_plan`, `plan_review`) stays on the login node so pairing-detection conflicts surface interactively before cluster jobs are dispatched.

### One-time setup

```bash
# Imperial RDS (PBS Pro):
export PMA_PBS_QUEUE=v1_throughput72
export PMA_PBS_PROJECT=<your project code>
export PMA_NOTIFY_EMAIL=<you@example.com>

# Generic SLURM site:
export PMA_SLURM_PARTITION=cpu
export PMA_SLURM_ACCOUNT=<your account>

# Optional — scale per-rule mem and walltime (default 1):
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

### Phase B — submit the unattended head-job (S1a → S3 + post_qc_review)

```bash
processing-muagent submit --config $CFG --executor pbs \
    --auto-approve --auto-approve-except post_qc_review \
    --auto-approve-except s7_clustering
```

The head-job stops at `post_qc_review` propose and emails `$PMA_NOTIFY_EMAIL` (if set). Poll progress:

```bash
processing-muagent status --watch --config $CFG
```

### Phase C — QC review + resume (S4 → S7 propose)

Review `deliverables/checkpoint/qc_review/qc_summary.md` and figures in `deliverables/checkpoint/qc_review/`. On paired runs, confirm the S3 doublet policy section. Revise thresholds if needed, then approve:

```bash
# Optionally revise, e.g.:
processing-muagent revise s2_atac_qc s2_atac_qc.tss_enrichment_min=1.5 --config $CFG

processing-muagent approve post_qc_review --config $CFG
processing-muagent submit --config $CFG --executor pbs \
    --auto-approve --auto-approve-except s7_clustering
```

The head-job stops again at `s7_clustering` propose.

### Phase D — resolution review + finish (S7 execute → S8 → manifest)

Open in any browser:

```
<run_dir>/deliverables/checkpoint/resolution_review/resolution_review.html
```

The accompanying `resolution_review.ipynb` is for power users who want to re-cluster at custom resolutions. Approve or revise:

```bash
processing-muagent approve s7_clustering --config $CFG
# OR:
processing-muagent revise s7_clustering s7_clustering.rna.resolution=1.2 --config $CFG

processing-muagent submit --config $CFG --executor pbs
```

### Foreground cluster mode (alternative to `submit`)

Keep Snakemake attached in your tmux session (lowest-latency approvals, no head-job queue time):

```bash
processing-muagent run --config $CFG --executor pbs
```

Snakemake dispatches per-rule cluster jobs and exits cleanly at unapproved gates. Re-invoke after each approval.

### Per-stage resources

Edit `workflow/resources.smk` to override mem/walltime/cpus. The table is the single source of truth for both PBS and SLURM profiles. OOM-killed jobs are retried once at 2× memory (`restart-times: 1`).

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

## Environment

Implementation developed against the `cell_annotation` micromamba env with pip-installed `muon`, `scrublet`, `leidenalg`, `snakemake`, `mudata`. The `workflow/envs/*.yaml` files are the canonical production conda definitions.

**Ambient-correction R dependency (optional).** S1a calls DecontX (`celda`) or SoupX (`SoupX`) via `Rscript`. If R / the requested package is not installed, S1a degrades to pass-through and records `s1a_ambient.method = "skipped_no_r"` in `parameters.yaml`. To enable:

```bash
Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("celda","SoupX"))'
```

SnapATAC2 function names (`pp.import_fragments`, `metrics.tsse`, `pp.add_tile_matrix`, `pp.select_features`, `tl.spectral`, `tl.leiden`, `tl.umap`) were selected for SnapATAC2 ≥ 2.6; verify against the installed version at execute time.
