# Processing-MuAgent

Multiome (scRNA-seq + scATAC-seq) preprocessing subagent. Takes filtered or raw 10x Genomics multiome outputs and performs QC, PCA (RNA) + neighbor graphs, clustering, and UMAP **per modality**, then **stops** before integration.

Supported workflow branches: `paired`, `separate`, `rna_only`, `atac_only`. Declare the branch up front with `Processing-MuAgent declare-branch`.

**Harness layout.** Agent manifest + identity: [`AGENT.md`](AGENT.md). Conversational procedures (progressive disclosure): [`agent/skills/`](agent/skills/) (start at `index.md`); policy + hard rules: [`agent/system_prompt.md`](agent/system_prompt.md); per-command tool contracts: [`agent/tools.md`](agent/tools.md). Cross-agent contracts (finding codes, state model, handoff schemas): [`../contracts/`](../contracts/). QC default values live once in `executor/defaults.py`. Repo overview: [`../README.md`](../README.md).

## Pipeline overview

**Stage order** (Snakemake DAG — `s0_ingest` is a single planning-compute job that loads the data once, validates it, assembles the preprocessing plan, and runs the QC threshold exploration; it emits `validation_report.json`, `preprocessing_plan.json`, and `qc_explore.json` that `plan_review` consumes):

```
  P1 context extraction → S0 ingest (load + validate + assemble plan + QC explore) → (CHECKPOINT 1) plan_review
  → S1a ambient RNA correction → S1 RNA QC → S2 ATAC QC → S3 doublets → (CHECKPOINT 2) post_qc_review
  → qc_handoff (post-QC Integration bundle; after QC approval, orthogonal to S4–S8)
  → S4 RNA normalization + HVG → S5 ATAC spectral embedding → S6 PCA (RNA) + neighbor graph
  → S7 clustering (fixed resolutions) → S8 UMAP → manifest
```

### User checkpoints (2)

Two deliberate pauses where you review deliverables and decide before heavy downstream work continues. Everything else — including clustering and UMAP — runs automatically once upstream artifacts exist and the relevant checkpoint is approved.


| #     | CLI name              | Internal stage   | When                   | What you decide                                                                                                                                                                                                                                                            |
| ----- | --------------------- | ---------------- | ---------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **1** | **Plan review**       | `plan_review`    | After S0, before S1    | Approve the preprocessing plan (`plan/plan_review_<run>.md`)                                                                                                                                                                                                             |
| **2** | **QC review**         | `post_qc_review` | After S3, before S4/S5 | Inspect QC figures in `deliverables/figures/` + `deliverables/qc/qc_review_<run>.md` (or `qc_summary_<run>.html`); revise **RNA/ATAC quality-filter thresholds** (or skip individual metrics entirely) and re-run if needed; on **paired** multiome, confirm the **union doublet removal policy** |

After QC approval the pipeline runs straight through to the final outputs: Leiden clustering uses **fixed per-modality resolutions (RNA = 0.7, ATAC = 0.5)**. To use different values, revise clustering resolutions at plan review, or manually re-run clustering on the final outputs using your chosen parameters.


## Workflow stages

### Planning (pre-QC)

- **P1 Context extraction** — Biological Context Report (organism, tissue, assay, DOIs) plus DOI-based prior-analysis extraction.
- **S0 Ingest** — Loads and validates the input files, determines the workflow branch, and prepares the materials for user review at **checkpoint #1**. It accepts both **filtered** and **raw** Cell Ranger matrices, automatically detecting RNA and ATAC formats (see tables below). For **paired** multiome runs, it checks if RNA and ATAC modalities share cell barcodes; if not, it switches to the `separate` branch and records the reason. S0 also performs a **QC threshold preview** — exploring data distributions, estimating per-cell quality cutoffs, reporting how many cells would be removed by each metric, and generating diagnostic histograms. Based on these assessments, it produces a preprocessing plan for user review before QC filtering begins.
  **Supported RNA input formats (`rna_path`):**

  | Format tag  | File pattern           | Notes                                                                                                                                                                                 |
  | ----------- | ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
  | `10x_h5`    | `*.h5`                 | Cell Ranger (ARC) HDF5; GEX features filtered automatically                                                                                                                           |
  | `10x_mex`   | directory              | 10x MEX bundle with `matrix.mtx[.gz]` + `barcodes.tsv[.gz]`                                                                                                                           |
  | `h5ad`      | `*.h5ad`               | AnnData; `.X` must contain raw integer counts                                                                                                                                         |
  | `dense_txt` | `*.txt.gz`, `*.tsv.gz` | Dense genes × cells tab-delimited matrix (common GEO supplementary format). Row 0 = cell-barcode header; rows 1+ = gene symbol + counts. Loaded in 500-gene chunks to bound peak RAM. |

  **Supported ATAC input formats (`atac_fragments_path`):**

  | Format tag              | File pattern                | Notes                                                                                                                                                                                                                                                                                                                                                 |
  | ----------------------- | --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
  | `fragments_tsv`         | `*.tsv.gz` + `*.tsv.gz.tbi` | Standard 5-column bgzipped fragments file (`chrom start end barcode count`); tabix index must be present                                                                                                                                                                                                                                              |
  | `bed4` *(auto-convert)* | `*.bed.gz`                  | 4-column BED (`chrom start end barcode`). S0 auto-converts to a standard 5-column `fragments.tsv.gz` using `zcat → awk → sort → bgzip → tabix`. The source file is **never modified**; the derived `.tsv.gz` + `.tbi` are written alongside it. Windows `\r\n` line endings are handled automatically. Requires `bgzip` and `tabix` (htslib) on PATH. |

- **plan_review** — Generates a summary at `deliverables/plan/plan_review_<run>.md` (including the QC threshold-preview tables + histograms) for the user to review. The workflow pauses here until approval, before any S1–S8 execute rule runs.

### Preprocessing

- **S1a Ambient RNA correction** — Default `method=auto` on RNA branches (SoupX if raw+filtered exist, else DecontX). Omitted on `atac_only`. Whether to run is confirmed by user at **plan review (#1)** depending on the study goal, inputs, and sample context (see [10x ambient RNA guide](https://www.10xgenomics.com/analysis-guides/introduction-to-ambient-rna-correction)).
- **S1 RNA QC** — The following per-cell thresholds are applied:
  - `total_counts` (total UMI counts per cell): Median Absolute Deviation (MAD)-derived threshold with lower floor = **500**
  - `n_genes_by_counts` (number of genes detected per cell): MAD-derived threshold with lower floor = **250**
  - `pct_counts_mt` (percentage of counts mapping to mitochondrial genes): MAD-derived upper threshold with lower floor = **5%** and ceiling = **20%**
  - `pct_counts_ribo` (percentage of counts mapping to ribosomal genes): upper ceiling = **50%**
- **S2 ATAC QC** — The following per-cell thresholds are applied:
  - `n_fragments` (number of fragments per cell): MAD-derived threshold with lower floor = **1,000**
  - `TSS_enrichment` (Transcription Start Site enrichment score): minimum = **1.5**, maximum = **50**
  - `nucleosome_signal` (nucleosome signal): default = **3**
  - `FRiP` (Fraction of Reads in Peaks): default = **0.2**

  Any MAD-derived bound above (`total_counts`, `n_genes_by_counts`, `pct_counts_mt`, `n_fragments`) can be pinned to an exact value with its `*_override` key; see **Flexible QC thresholds** below.

**Flexible QC thresholds**
Every RNA and ATAC QC metric can be **tightened/loosened**, **pinned to an exact value**, individually **skipped** (filter removed entirely), or **partially skipped** (upper or lower bound only removed) — at either **plan review** (checkpoint #1) or **QC review** (checkpoint #2). To pin a MAD-derived bound to a specific number, set its `*_override` key (e.g. `revise s1_rna_qc n_genes_min_override=300`); the MAD/floor derivation still runs and is drawn as a **grey** reference line while the chosen cutoff is drawn in **red** on the QC histograms. An override more permissive than the recommended floor/ceiling is still applied but flagged with a warning in the QC report.

- **S3 Doublets** — Per-modality doublet detection, then branch-specific reconciliation:
  - **RNA:** Scrublet (sparse-CSR input; `expected_doublet_rate ≈ 0.0008 × n_cells`, capped at 10%).
  - **RNA / ATAC:** fixed doublet score thresholds (defaults: RNA Scrublet 0.25, ATAC SnapATAC2 0.5; configurable via plan or `revise s3_doublets`).
  - **separate / single-modality branches:** Each modality is filtered independently by its own detector; per-modality calls are saved in `calls.parquet`.
  - **paired branch:** Also performs joint barcode alignment after doublet removal; the union doublet policy is confirmed at the **QC review checkpoint** (`deliverables/qc/qc_review_<run>.md`).
- **post_qc_review** — **QC review checkpoint (#2).** Propose-only gate between S3 and S6 PCA (RNA) + neighbor graph. Generates doublet histograms (with the chosen doublet threshold drawn as a red cutoff line), a cell-count waterfall (with counts labelled on bars), and `deliverables/qc/qc_review_<run>.md` — a plain-language summary of what each filter step did (MAD outlier bounds, MT/ribo ceilings, TSS enrichment, nucleosome signal, FRiP, union doublet policy). Each RNA/ATAC section opens with cells before filtering, retained, and removed. Revise quality-filter thresholds and re-run affected stages before approving. On approval, the large QC-only working files are automatically deleted to free storage (~2 GB/run): the QC matrices `rna_qc.h5ad`, `atac_qc.h5ad`, `atac_snap.h5ad`, the `qc_explore/atac_snap_explore.h5ad` import, the chr-normalised fragment caches `atac_fragments_cbf[_chrnorm].tsv.gz` (the single biggest artifact — reused across QC re-runs but dead once approved), and the S1a recompute caches (`tsne_coords_cache.parquet`, `cell_totals.parquet`). None is a declared Snakemake output or read by a post-gate stage, so deletion never triggers a re-run. Preserved: `qc_summary.json` markers, the QC-metrics parquets (the final S8 manifest reads `qc_metrics_post.parquet`), `rna_decontaminated.h5ad`, and all S3+ artifacts. (Per-stage scratch dirs — `_work_soupx`/`_work_decontx`, `macs3_tmp` — are removed by their own stage as soon as it finishes, not here.)
- **S4 RNA norm + HVG** — Log-normalize (`target_sum=1e4`) + HVG selection (`seurat_v3` on counts).
- **S5 ATAC spectral embedding and peak matrix export** — SnapATAC2 tile matrix (`bin_size=500`, unified with S3) → feature selection → spectral embedding. In parallel, exports a feature (cell-by-feature) matrix using this priority order for the peak coordinates:
  1. **User-supplied peaks** — `atac_peaks_path` in `run.yaml` → SnapATAC2 `make_peak_matrix` (`user_peaks` mode).
  2. **ARC peak matrix** — pre-called peaks from a combined Cell Ranger ARC `.h5` detected at S0 (`arc_h5` mode).
  3. **S2 pre-called peaks** — BED file written by S2 ATAC QC (MACS3 or ARC-derived) reused here; no redundant peak calling (`s2_peaks_macs3` / `s2_peaks_arc` mode).
  4. **Tile-matrix fallback** — verified SnapATAC2 tile matrix (`tile_matrix_fallback` mode), used only when no peak source is available.
  Spectral embedding in `obsm['X_spectral']` (with `X_lsi` as a backward-compat alias) is always computed from the tile matrix regardless of peak-export mode. When `drop_first=True`, the first component is removed before S6–S8.
- **S6 PCA (RNA) + neighbor graph** (`s6_neighbors`) — **RNA:** optional `sc.pp.scale`, then PCA; `n_pcs` from a chord-distance elbow on explained variance, capped at `rna_n_pcs_max`; nearest-neighbors on PCA space. **ATAC:** KNN graph on the S5 spectral embedding (`X_spectral` via `snap.pp.knn`). Artifact: `internal/artifacts/s6_neighbors/rna_neighbors.h5ad`.
- **S7 Clustering** — Leiden clustering at fixed per-modality resolutions (RNA = 0.7, ATAC = 0.5; `s7_clustering.rna_resolution` / `atac_resolution`). Runs automatically with no sweep and no checkpoint. Separate / single-modality branches: these become the final `leiden_rna` / `leiden_atac` labels. Paired branch: diagnostic per-modality labels for UMAP only (not joint embedding).
- **S8 UMAP** — Per-modality UMAP. **Paired** → `processed_<run>.h5mu`; **separate** → `rna_processed.h5ad` + `atac_processed.h5ad`. On the paired branch, S8 expects matching barcodes from S3; final assembly includes a defensive re-intersection logged only when it filters cells.
- **manifest** — `run_manifest.json` preprocessing handoff contract (v1.0.0), the review notebook, and `layout.json`.
- **qc_handoff** — After QC approval (gated on `post_qc_review.approved`; reads S3 post-doublet artifacts only). Writes `deliverables/qc/post_qc_<run>.h5mu` (post-QC, post-doublet, **un-normalized** cells for Integration-MuAgent; h5mu on all branches — paired, separate, rna_only, atac_only) and `deliverables/qc/post_qc_manifest.json` (schema `muagene.post_qc_handoff/1`). Orthogonal to S4–S8; `rule all` requires both this bundle and `run_manifest.json`. Independently buildable via `run --target qc_handoff`.

## Paired multiome

The paired branch admits three input shapes:

1. A single Cell Ranger ARC `.h5` (combined GEX + Peaks; barcodes match by construction).
2. Cell Ranger GEX `.h5` + ATAC fragments where cell barcodes match directly (or differ only by a `-N` / `_LIBRARY` suffix).
3. Independent GEX + ATAC pipelines whose barcodes live in different 10x whitelists — requires a 2-column TSV at `barcode_translation_path` (or `cell_metadata_path` with `rna_barcode` + `atac_barcode` columns) so S0 can rewrite ATAC barcodes into RNA space before QC.

In all cases, the final `processed_<run>.h5mu` contains only cells passing both RNA and ATAC QC with matching barcodes.

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

Doublets are detected separately in RNA and ATAC; their results are combined during S3 reconciliation using **union** policy: a cell is removed if either detector (RNA or ATAC) flags it as a doublet. Detectors are prone to false negatives, so union minimises contamination. Detector scores and flags are saved in `calls.parquet` for later review.

### Diagnostic vs final clustering

At S7, RNA-only (`leiden_rna`) and ATAC-only (`leiden_atac`) clustering are run separately for diagnostic comparison—not for joint clustering. In paired mode, both clusterings use the same intersected cell set from S3; in separate mode, they use their respective modality’s cells. Joint multimodal clustering (e.g., WNN or MOFA+) is not performed in this preprocessing workflow.

## Repository layout

```
Processing-MuAgent/
├── agent/               # system_prompt.md + skills/ (procedures) + tools.md (contracts)
├── config/              # example run configurations
├── executor/            # Python implementation (stages, methods, CLI, helpers)
│   ├── stages/          # per-stage scripts S0..S8 + post_qc_review + qc_handoff
│   ├── methods/         # MAD thresholds, doublet policy
│   ├── hpc.py           # site.config/hpc.env writers; submission delegated to Execution-MuAgent
│   └── specs.py         # stage metadata authoring (writes internal/stage_meta/)
├── workflow/            # Snakemake orchestration
│   ├── Snakefile        # localrules for planning + propose + manifest
│   ├── resources.smk    # per-stage mem/runtime/cpus + PROGRESS_TIMEOUT_HINT
│   ├── rules/           # per-stage propose/execute rule pairs + manifest
│   ├── envs/            # CPU env (processing.yaml -> muagene) + GPU container (muagene-gpu.def)
│   └── profiles/
│       └── slurm/       # SLURM snakemake profile
└── scripts/             # launch_runner.sh + head-job templates
```

## Running on HPC (SLURM)

On a cluster, heavy compute stages run as scheduler jobs (SLURM). The agent drives the workflow through the same checkpoints as local mode. Platform settings are gathered once per run via `configure-execution` and stored in `site.config` (the single source of truth); `hpc.env` is generated from it automatically. Everything else — init, submit, approve, revise — is handled via the CLI or the chat agent.

**Execution-MuAgent is a hard dependency for cluster submission.** `Processing-MuAgent submit` delegates rendering, submission, and monitoring to the sibling `Execution-MuAgent` package. If it is absent the command fails loudly. The **Install** bootstrap above (`init-machine`) installs it into the `muagene` env for you; on a machine set up some other way, `pip install -e Execution-MuAgent/` into the same env.

### One-time setup per run

```bash
# Probe available queues/partitions and suggest settings:
Processing-MuAgent hpc-info

# Write site.config (source of truth) and derived hpc.env.
# --confirmed-by-user records that the user approved this execution mode; run/submit
# refuse to launch any compute until it is set (applies to local mode too).
# Local:
Processing-MuAgent configure-execution --config $CFG --mode local --confirmed-by-user

# SLURM:
Processing-MuAgent configure-execution --config $CFG --mode slurm \
  --slurm-partition <partition> --slurm-account <account> --confirmed-by-user

# SLURM on GPU (cluster-only, for integration subagent — preprocessing is CPU-only):
Processing-MuAgent configure-execution --config $CFG --mode slurm \
  --slurm-partition <cpu-partition> --slurm-account <account> \
  --device gpu --gpu-partition gpu --gpu-gres gpu:A5000:1 \
  --confirmed-by-user
# SLURM --device gpu REQUIRES --gpu-image-uri (the GPU env is a container PULLED from that pinned
#  reference). It + --singularity-module + env manager auto-fill from ~/.muagene/machine.config if
#  you ran init-machine — otherwise pass --gpu-image-uri docker://<registry>/muagene-gpu:<tag>.
# Optional: --scratch <path> binds an extra node-local/fast path into the GPU container
#  (exported as PMA_GPU_BIND). The run directory and repo root are always bound.
```

GPU cluster infrastructure (container pull, partition/gres routing, bind contract) is
in place for the **integration subagent** (future). **Processing-MuAgent preprocessing
is CPU-only** — `_GPU_CAPABLE` in `workflow/resources.smk` is empty.

This writes:

- `deliverables/plan/config/site.config` — YAML platform description (consumed by Execution-MuAgent)
- `deliverables/plan/config/hpc.env` — shell snippet generated from site.config; source before `submit`

For larger datasets increase `--resources-scale` (e.g. `2` for ~30k cells, `4` for ~100k). Per-stage CPU, memory, and walltime defaults live in `workflow/resources.smk`; OOM-killed jobs are retried once at double memory.

### Requirements

- **SLURM:** Uses Snakemake's generic cluster executor plus `sbatch`/`sacct`/`squeue`; every heavy rule is submitted as a separate `sbatch` job, even when the Snakemake head-job itself runs under SLURM. Child jobscripts are sanitized before submission so Snakemake does not re-enable `storage-local-copies` on shared NFS (a common cause of jobs finishing Python quickly then hanging at `Storing output in storage.`).
- **Clean stage exit:** after writing its outputs, every cluster `<stage>_execute` rule calls `executor.cluster_exit.finalize_cluster_exit()` (`gc.collect()` + `os._exit(0)` when `SLURM_JOB_ID` is set). This terminates lingering h5py/HDF5 background threads that would otherwise keep the child process alive after the output is complete — the previous cause of jobs stuck RUNNING with `slurmstepd: error: Pid still in cpuset cgroup`. It is a no-op in local mode (job-id env vars unset).

### How the HPC run proceeds


| Step                                | Stages                      | Executes on                     | You                     |
| ----------------------------------- | --------------------------- | ------------------------------- | ----------------------- |
| Context                             | P1                          | Login node (cheap, interactive) | Fill biological context |
| load + validate + plan + QC explore | S0                          | Cluster                         | —                       |
| Checkpoint **#1**                   | plan_review                 | Login node                      | Review plan             |
| QC                                  | S1a → S1 → S2 → S3          | Cluster                         | —                       |
| Checkpoint **#2**                   | post_qc_review              | —                               | Review QC               |
| Integration handoff                 | qc_handoff                   | Cluster (SLURM) job; runs at QC approval via `submit --target qc_handoff` | — |
| Finish                              | S4 → S5 → S6 → S7 → S8 → manifest | Cluster | —                       |


**Execution mode must be user-confirmed before any compute runs.** Both `run` and `submit` hard-refuse to launch until `execution.user_confirmed=true` is recorded (via `configure-execution ... --confirmed-by-user`). This is a one-time gate enforced on fresh runs and resume sessions alike — the agent must confirm local vs HPC with the user and never auto-default. `run` additionally refuses when the mode is `slurm` (use `submit`). Once confirmed, the pipeline proceeds automatically.

**S0 execution mode:** in HPC mode (`execution.mode = slurm`), S0 is **always** submitted through Execution-MuAgent as a supervised cluster job (never run on the login node — its QC exploration needs 100+ GB). `submit` with no `--target` infers `plan_review_propose` as the planning target, which pulls P1 → S0 as dependencies and arms the plan-review gate in one head-job, before checkpoint #1; the supervision daemon monitors it, and you report its status with one-shot `hpc-status`. In local mode, run `run --target plan_review_propose` on this machine.

Each heavy stage runs as its own scheduler job, and only the two checkpoints above need `approve`. After `submit`, a background monitor is the sole watcher of your job; Processing-MuAgent follows **report-and-repoll** — it reports one-shot `hpc-status`, then re-polls on a non-blocking scheduled wakeup (~295s, the cadence printed on the `Next check:` line) and re-reports only when the state changes, until the job finishes or a review gate arms. You never have to ask for status by hand.

### Submit workflow

Source `deliverables/plan/config/hpc.env`, then use `Processing-MuAgent submit` (not `run`) to dispatch the Snakemake head-job. `**submit` auto-infers the Snakemake target** from run state — you do not need to pick `plan_review_propose`, `post_qc_review_propose`, or `all` manually. After each approval, run `submit` again and it stops at the next gate:


| Run state                                 | Inferred target          | Runs through                                                              |
| ----------------------------------------- | ------------------------ | ------------------------------------------------------------------------- |
| `plan_review` not approved                | `plan_review_propose`    | Fresh run: P1 → S0 → plan assembly + QC explore, then arms gate. Resume after S0: only the cheap propose rule runs. |
| `post_qc_review` not approved             | `post_qc_review_propose` | S1a → S3 + QC summary, then pauses                                        |
| `post_qc_review` approved                 | `all`                    | S4 → S6 → S7 clustering → S8 → manifest → final results (no further pause; `qc_handoff` already ran at QC approval and is skipped by Snakemake) |


Override with `--target <name>` only when debugging.

If a previous run was cancelled or killed, stale Snakemake locks under `<run_dir>/internal/snakemake/` can make the next submit fail immediately with `LockException`. Recover with:

```bash
Processing-MuAgent unlock --config $CFG
# or on submit:
Processing-MuAgent submit --config $CFG --executor slurm --unlock-stale-locks
```

```bash
CFG=<run_dir>/deliverables/plan/config/run.yaml

# Configure HPC settings (writes site.config + derived hpc.env). --confirmed-by-user
# records the user's approval of the mode; submit refuses to launch without it:
Processing-MuAgent configure-execution --config $CFG --mode slurm \
  --slurm-partition cpu-medium --slurm-account mylab --confirmed-by-user

source <run_dir>/deliverables/plan/config/hpc.env

# Submit the planning head-job via Execution-MuAgent (P1 → S0 → gate-arming);
# `submit` with no --target infers `plan_review_propose`:
Processing-MuAgent submit --config $CFG --executor slurm
# When it completes, render the plan-review deliverable and approve:
Processing-MuAgent plan-review --config $CFG  # also writes internal/stage_meta/
Processing-MuAgent approve plan_review --config $CFG

# First heavy QC batch (stops at QC review):
Processing-MuAgent submit --config $CFG --executor slurm

# After QC review — approve, then run qc_handoff immediately to create integration artifacts:
Processing-MuAgent approve qc_review --config $CFG
Processing-MuAgent submit --config $CFG --executor slurm --target qc_handoff
# When ready to run clustering + UMAP through to final outputs (no further pause):
Processing-MuAgent submit --config $CFG --executor slurm
```

After `submit` returns, the background monitor keeps watching your job — see **How cluster jobs run and stay supervised** below. For an unattended batch, pre-seed both checkpoints with `--auto-approve` on `run` or `submit`; to keep specific gates interactive, add `--auto-approve-except qc_review` (repeatable; accepts checkpoint aliases).

### How cluster jobs run and stay supervised

`Processing-MuAgent submit` hands your job to the sibling **Execution-MuAgent** package, which prepares the scheduler script, submits the job, and then runs a **background monitor** that watches it for its whole lifetime. Install Execution-MuAgent first — `submit` fails loudly without it. `submit` returns in about 90 seconds, once the job is confirmed accepted; the monitor keeps running after that and is the only thing watching your job.

**The monitor is a safety layer, not optional.** It watches for stalled or hung jobs and cancels one only when the evidence is conclusive — so a genuinely stuck job doesn't sit burning cluster time, while a slow-but-alive job is left alone. It also cleans up the opposite case: a job whose work has finished but whose process lingers is cancelled so it can't burn its allocation to walltime. As each stage finishes it also confirms the stage actually produced a complete, readable output. If the monitor dies, your cluster job keeps running but is no longer protected. Check the monitor's liveness and your job's health at any time with one-shot `hpc-status` (no poll loop):

```bash
Processing-MuAgent hpc-status --config $CFG
```

**Keeping the monitor alive across logout.** It runs in the background and survives SSH disconnects on most systems. Some clusters kill all your processes when you log out; on those, start your session inside `tmux` or `screen` before `submit`. If the monitor does die while the job is still running, restart it without resubmitting:

```bash
Processing-MuAgent supervisor-restart --config $CFG
```

**When a job fails or is rejected.** Execution-MuAgent never contacts you directly and never resubmits — it records what it found and stops. `hpc-status` shows the result: if the scheduler rejected the job (usually a wrong partition, account, or walltime), re-run `configure-execution` with corrected settings; for any other failure, fix the cause and re-run `submit`. Cluster jobs only ever run through `submit` — there is no manual submission path.

Each `submit` first **archives the previous run's Snakemake logs** to `<run_dir>/internal/snakemake/.snakemake/archive/run_<timestamp>/` (a move, so history is preserved). This keeps the live log dirs scoped to the current run, so `hpc-status` reads stage state from the new job's logs only and never reports a phantom failure from the prior run during the new job's PENDING window. The root-cause line for the failure that prompted the resubmit is in the archived logs (and was already surfaced via `hpc-status` before you resubmitted).

For how the monitor works internally — how it decides a job is stalled, how it verifies outputs, and the diagnostics it records — see `**Execution-MuAgent/README.md`**, the package that owns it.

## Run directory layout

After `init`, only `deliverables/plan/` exists under deliverables. The `figures/`, `qc/`, and `results/` folders are created when the pipeline first writes into them.

```
<run_dir>/
  deliverables/
    plan/
      config/
        run.yaml                  ← canonical config (use this for all CLI calls)
        biological_context.md     ← Biological Context Report
        hpc.env                   ← source before submit/run on cluster
        site.config               ← YAML platform description (consumed by Execution-MuAgent)
      context_summary.md          ← P1 output
      plan_review_<run>.md        ← plan review gate (summary + parameter appendix)
      plan_summary_<run>.html     ← download-friendly web version of plan_review (figures embedded)
    figures/                      ← all pipeline figures (created at first plot)
    qc/                           ← QC checkpoint (#2): reports + post-QC Integration handoff
      qc_review_<run>.md          ← before approval
      qc_summary_<run>.html
      post_qc_<run>.h5mu          ← after approval (qc_handoff; all branches)
      post_qc_manifest.json
    results/                      ← final deliverables (S8 + manifest)
      processed_<run>.h5mu         ← or rna/atac_processed.h5ad (separate branch)
      review_processed_<run>.{ipynb,py}
      run_manifest.json            ← manifest (preprocessing handoff)
      layout.json
  internal/
    artifacts/sN_<stage>/         ← intermediate stage outputs
    stage_meta/<stage>.yaml       ← per-stage metadata (resources, I/O, timeout hint) — not a submission contract
    stage_meta/head_job.yaml      ← head-job submission spec (written by submit, read by Execution-MuAgent)
    proposals/                    ← optional <stage>.yaml (mainly checkpoint review artifacts)
    checkpoints/                  ← plan_review, post_qc_review .approved only
    hpc_monitor/
      submissions.jsonl           ← append-only job registration log
      latest_submission.json      ← most recent submission record (used by supervisor-restart)
      execution_manifest.jsonl    ← per-submit record (stage, job_id, script, outputs)
      latest_report.md            ← Execution-MuAgent investigation reports and confirmed-dead verdicts
      latest_snapshot.json        ← full snapshot + monitor_state (health, silence, investigation)
      monitor.pid                 ← PID of the running supervision daemon (removed when daemon exits)
      monitor.log                 ← symlink to the latest monitor_<timestamp>.log
      monitor_<timestamp>.log     ← daemon output for each submit or supervisor-restart
    parameters.yaml
    state.yaml
    log.jsonl
```

## CLI

Commands below use `$CFG` = `<run_dir>/deliverables/plan/config/run.yaml` (written by `init`).

### Install

MuAgene is set up on a fresh machine by **Execution-MuAgent** (it owns infrastructure). Clone both repos as siblings, then run **one** bootstrap command — it creates the single integrated `muagene` env (science stack + both agent CLIs) from the committed conda-lock lock and installs both packages into it. Do **not** create the conda env by hand.

```bash
# With Processing-MuAgent and Execution-MuAgent cloned as siblings:
bash /path/to/Execution-MuAgent/scripts/bootstrap.sh --processing-repo /path/to/Processing-MuAgent
conda activate muagene
```

For GPU (`--device both --gpu-image-uri docker://<registry>/muagene-gpu:<tag>`) and the `~/.muagene/machine.config` profile, see `Execution-MuAgent/README.md`.

### Configure and scaffold a run

Edit `config/run.example.yaml` (at minimum `run_dir`, `genome_assembly`, and modality paths). Optional: `rna_raw_path`, `atac_peaks_path`, `barcode_translation_path`, `cell_metadata_path`, `biological_context_path`.

```bash
Processing-MuAgent init --config config/run.example.yaml
CFG=<run_dir>/deliverables/plan/config/run.yaml
Processing-MuAgent declare-branch paired --config $CFG   # paired | separate | rna_only | atac_only
```

`init` creates `<run_dir>/`, copies config to `deliverables/plan/config/run.yaml`, and writes the Biological Context Report template at `deliverables/plan/config/biological_context.md`.

### Command reference


| Command               | Purpose                                                                                                                                                                                        |
| --------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `init`                | Create run directory scaffold                                                                                                                                                                  |
| `declare-branch`      | Record workflow branch in `parameters.yaml`                                                                                                                                                    |
| `configure-execution` | Set `execution.mode` + `execution.user_confirmed` (`--confirmed-by-user`); write `site.config` (platform source of truth) and derived `hpc.env`. `run`/`submit` refuse compute until confirmed |
| `hpc-info`            | Probe SLURM queues, partitions, accounts on the login node                                                                                                                                 |
| `run`                 | Foreground Snakemake, **local-only** (local-mode execution + login-node localrules). No cluster path — use `submit` for SLURM                                                              |
| `submit`              | **The only cluster-execution path.** Submit head-job via Execution-MuAgent (hard dependency); starts background supervision daemon; infers phase target                                        |
| `supervisor-restart`  | Restart the supervision daemon without resubmitting — use when the daemon died mid-run (SSH drop, site reboot, OOM) but the cluster job is still active                                        |
| `status`              | Per-step pipeline state (S1a–S8 + review gates); `--watch` polls until actionable                                                                                                              |
| `hpc-status`          | One-shot report of job health, supervisor liveness, monitor findings, and per-step state (no poll loop); warns if supervision is offline while the cluster job is still running                |
| `approve`             | Write `internal/checkpoints/<stage>.approved` (human checkpoints only)                                                                                                                         |
| `plan-review`         | Render `plan_review.md`; also writes per-stage metadata to `internal/stage_meta/`                                                                                                              |
| `revise`              | Update one or more parameters in `parameters.yaml` and reset a checkpoint to awaiting. Used to tune, tighten, loosen, **pin to an exact value** (`*_override`), or **skip** individual QC metrics — see **Flexible QC thresholds**       |
| `regenerate-locks`    | Regenerate the CPU conda-lock lockfile from `workflow/envs/processing.yaml` after editing dependencies (needs `pip install '.[dev]'`; commit the refreshed lock)                                |
| `unlock`              | Remove stale Snakemake locks after a cancelled/killed run                                                                                                                                      |
| `propose`             | Run a single `*_propose` rule (optional; not required for the main pipeline)                                                                                                                   |


**Approve aliases:** `qc_review` → `post_qc_review`. `revise` accepts a short `<param>=<value>` form — the stage prefix is auto-added (e.g. `revise s7_clustering rna_resolution=1.2` stores `s7_clustering.rna_resolution`). The full `<stage>.<param>=<value>` form is also accepted.

### Local workflow

Planning and QC stages run automatically. Snakemake stops only at the two checkpoints; after QC approval it runs straight through clustering and UMAP to the final outputs.

```bash
CFG=<run_dir>/deliverables/plan/config/run.yaml

# Confirm execution mode once (required before any compute — run/submit refuse otherwise):
Processing-MuAgent configure-execution --config $CFG --mode local --confirmed-by-user

# Option A — pause at each checkpoint (recommended first time):
Processing-MuAgent run --config $CFG
Processing-MuAgent approve plan_review --config $CFG
Processing-MuAgent run --config $CFG
Processing-MuAgent approve qc_review --config $CFG
Processing-MuAgent run --config $CFG --target qc_handoff   # integration artifacts; runs immediately at approval
# When ready to proceed with S4 → S8 → manifest (no further pause):
Processing-MuAgent run --config $CFG

# Option B — pre-seed both checkpoints (unattended Snakemake; you still review outputs):
Processing-MuAgent run --config $CFG --auto-approve

# Option C — unattended except one gate (example: keep QC review interactive):
Processing-MuAgent run --config $CFG --auto-approve --auto-approve-except qc_review

Processing-MuAgent status --watch --config $CFG
```

`run` requires a filled Biological Context Report unless you pass `--no-context`. After `revise`, approve the affected checkpoint again before resuming.

### HPC workflow

After checkpoint **#1**, use `submit` instead of foreground `run`. See **Running on HPC → Submit workflow** for the resume loop, `unlock`, and `--unlock-stale-locks`.

```bash
source <run_dir>/deliverables/plan/config/hpc.env
Processing-MuAgent submit --config $CFG --executor slurm
# submit returns within ~90 s; the supervision daemon keeps running in the background as the sole monitor.
Processing-MuAgent hpc-status --config $CFG   # one-shot: job health, supervisor liveness, findings, per-step state, and a `Next check:` re-poll cadence (report-and-repoll)

# If the daemon dies mid-run (SSH drop, site reboot) but the cluster job is still running:
Processing-MuAgent supervisor-restart --config $CFG
```

### Optional debugging

```bash
Processing-MuAgent revise s7_clustering rna_resolution=1.2 --config $CFG
Processing-MuAgent run --config $CFG --no-context
Processing-MuAgent hpc-info
Processing-MuAgent propose post_qc_review --config $CFG
```

## Environment

Environment setup/management is owned by **Execution-MuAgent** (it owns the non-scientific runtime layer). This repo only *authors* the definitions in `workflow/envs/`: `processing.yaml` (CPU source-of-truth → the committed conda-lock lock `processing.linux-64.lock`) and `muagene-gpu.def` (the GPU container recipe). The per-device provider + paths live in one committed file, `workflow/envs/manifest.yaml`, read by both agents; `site.config`'s `environments:` section is generated from it. There is no manual `conda env create` step — the env is provisioned from the lock.

On a fresh machine, `Execution-MuAgent init-machine` does the whole setup (see **Install** above and `Execution-MuAgent/README.md`). Then:

- **CPU env** = the conda-lock lock — **linux-only** (a non-linux host fails loud with `platform_unsupported`, never a silent solve). Edited `processing.yaml`? Regenerate the lock and commit it: `Processing-MuAgent regenerate-locks` (needs `pip install '.[dev]'`). `validate-env`/`submit` fail loud (`lock_stale_vs_yaml`) when the YAML is newer than the lock.
- **GPU env** = a pinned container image **pulled** from a registry — built + published centrally from `muagene-gpu.def` (see `scripts/build_and_push_gpu_image.sh`), **never built on a target machine**.
- `submit` auto-provisions a missing/stale env before launching (policy=auto). GPU env is for the integration subagent (future); preprocessing is CPU-only.

Per-run provisioning/validation is optional once the machine is bootstrapped (`init-machine` already did it); to (re)provision for a specific run's `site.config`:

```bash
Execution-MuAgent provision-env --site-config <run>/deliverables/plan/config/site.config --repo-root . --device both
Execution-MuAgent validate-env  --site-config <run>/deliverables/plan/config/site.config --repo-root .
```

