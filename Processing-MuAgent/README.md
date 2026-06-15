# Processing-MuAgent

Multiome (scRNA-seq + scATAC-seq) preprocessing subagent. Takes filtered or raw 10x Genomics multiome outputs and performs QC, PCA (RNA) + neighbor graphs, clustering, and UMAP **per modality**, then **stops** before integration.

Supported workflow branches: `paired`, `separate`, `rna_only`, `atac_only`. Declare the branch up front with `Processing-MuAgent declare-branch`.

## Pipeline overview

**Stage order** (Snakemake DAG — `s0_ingest` is a single planning-compute job that loads the data once, validates it, assembles the preprocessing plan, and runs the QC threshold exploration; it emits `validation_report.json`, `preprocessing_plan.json`, and `qc_explore.json` that `plan_review` consumes):

```
  P1 context extraction → S0 ingest (load + validate + assemble plan + QC explore) → (CHECKPOINT 1) plan_review
  → S1a ambient RNA correction → S1 RNA QC → S2 ATAC QC → S3 doublets → (CHECKPOINT 2) post_qc_review
  → S4 RNA normalization + HVG → S5 ATAC spectral embedding → S6 PCA (RNA) + neighbor graph
  → S7 clustering + (CHECKPOINT 3) resolution_review → S8 UMAP → outputs
```

### User checkpoints (3)

Three deliberate pauses where you review deliverables and decide before heavy downstream work continues. Everything else runs automatically once upstream artifacts exist and the relevant checkpoint is approved.


| #     | CLI name              | Internal stage   | When                   | What you decide                                                                                                                                                                                                                                                            |
| ----- | --------------------- | ---------------- | ---------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **1** | **Plan review**       | `plan_review`    | After S0, before S1    | Approve the preprocessing plan (`plan/summary/plan_review.md`)                                                                                                                                                                                                             |
| **2** | **QC review**         | `post_qc_review` | After S3, before S4/S5 | Inspect QC figures in `deliverables/figures/` + `checkpoints/qc_review/qc_review_<run>.md`; revise **RNA/ATAC quality-filter thresholds** (or skip individual metrics entirely) and re-run if needed; on **paired** multiome, confirm the **union doublet removal policy** |
| **3** | **Resolution review** | `s7_clustering`  | After S6, before S8    | Choose Leiden resolution per modality from sweep metrics (`checkpoints/resolution_review/`). **Separate / single-modality:** sets **final** cluster labels. **Paired:** **diagnostic** per-modality labels for UMAP only (not joint embedding)                             |


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

- **plan_review** — Generates a summary at `deliverables/plan/summary/plan_review.md` (including the QC threshold-preview tables + histograms) for the user to review. The workflow pauses here until approval, before any S1–S8 execute rule runs.

### Preprocessing

- **S1a Ambient RNA correction** — Default `method=auto` on RNA branches (SoupX if raw+filtered exist, else DecontX). Omitted on `atac_only`. Whether to run is confirmed by user at **plan review (#1)** depending on the study goal, inputs, and sample context (see [10x ambient RNA guide](https://www.10xgenomics.com/analysis-guides/introduction-to-ambient-rna-correction)).
- **S1 RNA QC** — The following per-cell thresholds are applied:
  - `total_counts` (total UMI counts per cell): Median Absolute Deviation (MAD)-derived threshold with lower floor = **500**
  - `n_genes_by_counts` (number of genes detected per cell): MAD-derived threshold with lower floor = **250**
  - `pct_counts_mt` (percentage of counts mapping to mitochondrial genes): MAD-derived upper threshold with lower floor = **5%** and ceiling = **20%**
  - `pct_counts_ribo` (percentage of counts mapping to ribosomal genes): upper ceiling = **50%**
- **S2 ATAC QC** — The following per-cell thresholds are applied:
  - `n_fragments` (number of fragments per cell): MAD-derived threshold with lower floor = **1,500**
  - `TSS_enrichment` (Transcription Start Site enrichment score): minimum = **1.5**, maximum = **50**
  - `nucleosome_signal` (nucleosome signal): defult = **3**
  - `FRiP` (Fraction of Reads in Peaks): defult = **0.25**

**Flexible QC thresholds**
Every RNA and ATAC QC metric can be **tightened/loosened**, individually **skipped** (filter removed entirely), or **partially skipped** (upper or lower bound only removed) — at either **plan review** (checkpoint #1) or **QC review** (checkpoint #2).

- **S3 Doublets** — Per-modality doublet detection, then branch-specific reconciliation:
  - **RNA:** Scrublet (sparse-CSR input; `expected_doublet_rate ≈ 0.0008 × n_cells`, capped at 10%).
  - **RNA / ATAC:** fixed doublet score thresholds (defaults: RNA Scrublet 0.25, ATAC SnapATAC2 0.5; configurable via plan or `revise s3_doublets`).
  - **separate / single-modality branches:** Each modality is filtered independently by its own detector; per-modality calls are saved in `calls.parquet`.
  - **paired branch:** Also performs joint barcode alignment after doublet removal; the union doublet policy is confirmed at the **QC review checkpoint** (`checkpoints/qc_review/qc_review_<run>.md`).
- **post_qc_review** — **QC review checkpoint (#2).** Propose-only gate between S3 and S6 PCA (RNA) + neighbor graph. Generates doublet histograms, a cell-count waterfall (with counts labelled on bars), and `checkpoints/qc_review/qc_review_<run>.md` — a plain-language summary of what each filter step did (MAD outlier bounds, MT/ribo ceilings, TSS enrichment, nucleosome signal, FRiP, union doublet policy). Each RNA/ATAC section opens with cells before filtering, retained, and removed. Revise quality-filter thresholds and re-run affected stages before approving. On approval, the large intermediate QC objects `rna_qc.h5ad`, `atac_qc.h5ad`, and `atac_snap.h5ad` are automatically deleted to free storage; `qc_summary.json` files and QC metrics parquets are preserved for report generation and threshold revision.
- **S4 RNA norm + HVG** — Log-normalize (`target_sum=1e4`) + HVG selection (`seurat_v3` on counts).
- **S5 ATAC spectral embedding and peak matrix export** — SnapATAC2 tile matrix (`bin_size=500`, unified with S3) → feature selection → `snap.tl.spectral` (Laplacian eigenmaps with IDF feature weights; not classical TF-IDF + SVD LSI). In parallel, exports a feature (cell-by-feature) matrix using this priority order for the peak coordinates:
  1. **User-supplied peaks** — `atac_peaks_path` in `run.yaml` → SnapATAC2 `make_peak_matrix` (`user_peaks` mode).
  2. **ARC peak matrix** — pre-called peaks from a combined Cell Ranger ARC `.h5` detected at S0 (`arc_h5` mode).
  3. **S2 pre-called peaks** — BED file written by S2 ATAC QC (MACS3 or ARC-derived) reused here; no redundant peak calling (`s2_peaks_macs3` / `s2_peaks_arc` mode).
  4. **Tile-matrix fallback** — verified SnapATAC2 tile matrix (`tile_matrix_fallback` mode), used only when no peak source is available.
  Spectral embedding in `obsm['X_spectral']` (with `X_lsi` as a backward-compat alias) is always computed from the tile matrix regardless of peak-export mode. When `drop_first=True`, the first component is removed before S6–S8.
- **S6 PCA (RNA) + neighbor graph** (`s6_neighbors`) — **RNA:** optional `sc.pp.scale`, then PCA; `n_pcs` from a chord-distance elbow on explained variance, capped at `rna_n_pcs_max`; nearest-neighbors on PCA space. **ATAC:** KNN graph on the S5 spectral embedding (`X_spectral` via `snap.pp.knn`). Artifact: `internal/artifacts/s6_neighbors/rna_neighbors.h5ad`.
- **S7 Clustering** — Leiden resolution sweep with per-modality grid and stable-region knee picker. **Resolution review checkpoint (#3):** `checkpoints/resolution_review/resolution_review.html` / `.ipynb`. Separate branch: chosen resolutions become final labels. Paired branch: diagnostic per-modality labels for UMAP only.
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

Doublets are detected separately in RNA and ATAC; their results are combined during S3 reconciliation using **union** policy: a cell is removed if either detector (RNA or ATAC) flags it as a doublet. Detectors are prone to false negatives, so union minimises contamination. Detector scores and flags are saved in `calls.parquet` for later review.

### Diagnostic vs final clustering

At S7, RNA-only (`leiden_rna`) and ATAC-only (`leiden_atac`) clustering are run separately for diagnostic comparison—not for joint clustering. In paired mode, both clusterings use the same intersected cell set from S3; in separate mode, they use their respective modality’s cells. Joint multimodal clustering (e.g., WNN or MOFA+) is not performed in this preprocessing workflow.

## Repository layout

```
Processing-MuAgent/
├── agent/               # chat-runtime prompts (system_prompt, interaction_flow)
├── config/              # example run configurations
├── executor/            # Python implementation (stages, methods, CLI, helpers)
│   ├── stages/          # per-stage scripts S0..S8 + post_qc_review
│   ├── methods/         # MAD thresholds, resolution sweep, doublet policy
│   ├── hpc.py           # site.config/hpc.env writers; submission delegated to Execution-MuAgent
│   └── specs.py         # stage metadata authoring (writes internal/stage_meta/)
├── workflow/            # Snakemake orchestration
│   ├── Snakefile        # localrules for planning + propose + manifest
│   ├── resources.smk    # per-stage mem/runtime/cpus + PROGRESS_TIMEOUT_HINT
│   ├── rules/           # per-stage propose/execute rule pairs + manifest
│   ├── envs/            # conda env (mirrors `grn`)
│   └── profiles/
│       ├── pbs/         # PBS Pro snakemake profile
│       └── slurm/       # SLURM snakemake profile
└── scripts/             # launch_runner.sh + head-job templates
```

## Running on HPC (PBS Pro or SLURM)

On a cluster, heavy compute stages run as scheduler jobs (PBS Pro or SLURM). The agent drives the workflow through the same checkpoints as local mode. Platform settings are gathered once per run via `configure-execution` and stored in `site.config` (the single source of truth); `hpc.env` is generated from it automatically. Everything else — init, submit, approve, revise — is handled via the CLI or the chat agent.

**Execution-MuAgent is a hard dependency for cluster submission.** `Processing-MuAgent submit` delegates rendering, submission, and monitoring to the sibling `Execution-MuAgent` package. If it is absent the command fails loudly. Install it first: `pip install -e Execution-MuAgent/`.

### One-time setup per run

```bash
# Probe available queues/partitions and suggest settings:
Processing-MuAgent hpc-info

# Write site.config (source of truth) and derived hpc.env.
# --confirmed-by-user records that the user approved this execution mode; run/submit
# refuse to launch any compute until it is set (applies to local mode too).
# Local:
Processing-MuAgent configure-execution --config $CFG --mode local --confirmed-by-user

# PBS Pro:
Processing-MuAgent configure-execution --config $CFG --mode pbs \
  --pbs-queue <queue> --pbs-project <project> --confirmed-by-user

# SLURM:
Processing-MuAgent configure-execution --config $CFG --mode slurm \
  --slurm-partition <partition> --slurm-account <account> --confirmed-by-user
```

This writes:

- `deliverables/plan/config/site.config` — YAML platform description (consumed by Execution-MuAgent)
- `deliverables/plan/config/hpc.env` — shell snippet generated from site.config; source before `submit`

For larger datasets increase `--resources-scale` (e.g. `2` for ~30k cells, `4` for ~100k). Per-stage CPU, memory, and walltime defaults live in `workflow/resources.smk`; OOM-killed jobs are retried once at double memory.

### Requirements

- **SLURM:** Uses Snakemake's generic cluster executor plus `sbatch`/`sacct`/`squeue`; every heavy rule is submitted as a separate `sbatch` job, even when the Snakemake head-job itself runs under SLURM. Child jobscripts are sanitized before submission so Snakemake does not re-enable `storage-local-copies` on shared NFS (a common cause of jobs finishing Python quickly then hanging at `Storing output in storage.`).
- **Clean stage exit:** after writing its outputs, every cluster `<stage>_execute` rule calls `executor.cluster_exit.finalize_cluster_exit()` (`gc.collect()` + `os._exit(0)` when `SLURM_JOB_ID`/`PBS_JOBID` is set). This terminates lingering h5py/HDF5 background threads that would otherwise keep the child process alive after the output is complete — the previous cause of jobs stuck RUNNING with `slurmstepd: error: Pid still in cpuset cgroup`. It is a no-op in local mode (job-id env vars unset).

### How the HPC run proceeds


| Step                                | Stages                      | Executes on                     | You                     |
| ----------------------------------- | --------------------------- | ------------------------------- | ----------------------- |
| Context                             | P1                          | Login node (cheap, interactive) | Fill biological context |
| load + validate + plan + QC explore | S0                          | Cluster                         | —                       |
| Checkpoint **#1**                   | plan_review                 | Login node                      | Review plan             |
| QC                                  | S1a → S1 → S2 → S3          | Cluster                         | —                       |
| Checkpoint **#2**                   | post_qc_review              | —                               | Review QC               |
| PCA + neighbors + clustering        | S4 → S5 → S6 → S7 (sweep)   | Cluster                         | —                       |
| Checkpoint **#3**                   | s7_clustering               | —                               | Review resolution       |
| Finish                              | S7 (labels) → S8 → manifest | Cluster                         | —                       |


**Execution mode must be user-confirmed before any compute runs.** Both `run` and `submit` hard-refuse to launch until `execution.user_confirmed=true` is recorded (via `configure-execution ... --confirmed-by-user`). This is a one-time gate enforced on fresh runs and resume sessions alike — the agent must confirm local vs HPC with the user and never auto-default. `run` additionally refuses when the mode is `pbs`/`slurm` (use `submit`). Once confirmed, the pipeline proceeds automatically.

**S0 execution mode:** in HPC mode (`execution.mode = pbs | slurm`), S0 is **always** submitted through Execution-MuAgent as a supervised cluster job (never run on the login node — its QC exploration needs 100+ GB). `submit` with no `--target` infers `s0_ingest_execute` as the planning target and dispatches it as the first cluster job, before checkpoint #1; the supervision daemon monitors it, and you report its status with one-shot `hpc-status`. In local mode, S0 runs in the foreground on this machine via `run`.

Each heavy stage runs as its own scheduler job, and only the three checkpoints above need `approve`. After `submit`, a background monitor is the sole watcher of your job; Processing-MuAgent follows **report-and-repoll** — it reports one-shot `hpc-status`, then re-polls on a non-blocking scheduled wakeup (~295s, the cadence printed on the `Next check:` line) and re-reports only when the state changes, until the job finishes or a review gate arms. You never have to ask for status by hand.

### Submit workflow

Source `deliverables/plan/config/hpc.env`, then use `Processing-MuAgent submit` (not `run`) to dispatch the Snakemake head-job. `**submit` auto-infers the Snakemake target** from run state — you do not need to pick `s0_ingest_execute`, `post_qc_review_propose`, `s7_clustering_propose`, or `all` manually. After each approval, run `submit` again and it stops at the next gate:


| Run state                                 | Inferred target          | Runs through                                                              |
| ----------------------------------------- | ------------------------ | ------------------------------------------------------------------------- |
| planning not done                         | `s0_ingest_execute`      | load + validate + assemble plan + QC explore, then pauses for plan review |
| planning done, `plan_review` not approved | `plan_review_propose`    | renders the plan-review deliverable, then pauses                          |
| `post_qc_review` not approved             | `post_qc_review_propose` | S1a → S3 + QC summary, then pauses                                        |
| `s7_clustering` not approved              | `s7_clustering_propose`  | S4 → S6 + resolution sweep, then pauses                                   |
| All approved                              | `all`                    | S7 labels → S8 → manifest                                                 |


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

# Submit the planning job (load + validate + plan + QC explore) via Execution-MuAgent;
# `submit` with no --target infers `s0_ingest_execute` as the planning target:
Processing-MuAgent submit --config $CFG --executor slurm
# When it completes, render the plan-review deliverable and approve:
Processing-MuAgent plan-review --config $CFG  # also writes internal/stage_meta/
Processing-MuAgent approve plan_review --config $CFG

# First heavy QC batch (stops at QC review):
Processing-MuAgent submit --config $CFG --executor slurm

# After QC review:
Processing-MuAgent approve qc_review --config $CFG
Processing-MuAgent submit --config $CFG --executor slurm

# After resolution review:
Processing-MuAgent approve resolution_review --config $CFG
Processing-MuAgent submit --config $CFG --executor slurm
```

After `submit` returns, the background monitor keeps watching your job — see **How cluster jobs run and stay supervised** below. For an unattended batch, pre-seed all three checkpoints with `--auto-approve` on `run` or `submit`; to keep specific gates interactive, add `--auto-approve-except qc_review` (repeatable; accepts checkpoint aliases).

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

After `init`, only `deliverables/plan/` exists under deliverables. The `figures/`, `checkpoints/`, and `results/` folders are created when the pipeline first writes into them.

```
<run_dir>/
  deliverables/
    plan/
      config/
        run.yaml                  ← canonical config (use this for all CLI calls)
        biological_context.md     ← Biological Context Report
        hpc.env                   ← source before submit/run on cluster
        site.config               ← YAML platform description (consumed by Execution-MuAgent)
      summary/
        context_summary.md        ← P1 output
        plan_review_<run>.md      ← plan review gate (summary + parameter appendix)
        plan_summary_<run>.html   ← download-friendly web version of plan_review (figures embedded)
    figures/                      ← all pipeline figures (created at first plot)
    checkpoints/                  ← review reports (created at first checkpoint)
      qc_review/                    ← QC review (#2): qc_review_<run>.md + qc_summary_<run>.html
      resolution_review/          ← resolution_summary.md + resolution_review.{html,ipynb}
    results/                      ← final deliverables (created at S8/manifest; data + manifest)
      processed.h5mu              ← or rna/atac_processed.h5ad (separate branch)
      review_processed_h5mu.{ipynb,py}
      qc_summary.md               ← final QC summary
      run_manifest.json           ← handoff artifact
      layout.json
  internal/
    artifacts/sN_<stage>/         ← intermediate stage outputs
    stage_meta/<stage>.yaml       ← per-stage metadata (resources, I/O, timeout hint) — not a submission contract
    stage_meta/head_job.yaml      ← head-job submission spec (written by submit, read by Execution-MuAgent)
    proposals/                    ← optional <stage>.yaml (mainly checkpoint review artifacts)
    checkpoints/                  ← plan_review, post_qc_review, s7_clustering .approved only
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

```bash
cd /path/to/Processing-MuAgent
micromamba env create -n grn -f workflow/envs/processing.yaml   # once per machine
micromamba activate grn
pip install -e .
```

### Configure and scaffold a run

Edit `config/run.example.yaml` (at minimum `run_dir`, `genome_assembly`, `study_goal`, and modality paths). Optional: `rna_raw_path`, `atac_peaks_path`, `barcode_translation_path`, `cell_metadata_path`, `biological_context_path`.

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
| `hpc-info`            | Probe PBS/SLURM queues, partitions, accounts on the login node                                                                                                                                 |
| `run`                 | Foreground Snakemake, **local-only** (local-mode execution + login-node localrules). No cluster path — use `submit` for PBS/SLURM                                                              |
| `submit`              | **The only cluster-execution path.** Submit head-job via Execution-MuAgent (hard dependency); starts background supervision daemon; infers phase target                                        |
| `supervisor-restart`  | Restart the supervision daemon without resubmitting — use when the daemon died mid-run (SSH drop, site reboot, OOM) but the cluster job is still active                                        |
| `status`              | Per-step pipeline state (S1a–S8 + review gates); `--watch` polls until actionable                                                                                                              |
| `hpc-status`          | One-shot report of job health, supervisor liveness, monitor findings, and per-step state (no poll loop); warns if supervision is offline while the cluster job is still running                |
| `approve`             | Write `internal/checkpoints/<stage>.approved` (human checkpoints only)                                                                                                                         |
| `plan-review`         | Render `plan_review.md`; also writes per-stage metadata to `internal/stage_meta/`                                                                                                              |
| `revise`              | Update one or more parameters in `parameters.yaml` and reset a checkpoint to awaiting. Used to tune, tighten, loosen, or **skip** individual QC metrics — see **Flexible QC thresholds**       |
| `resolution-compare`  | Side-by-side resolution comparison figures (optional)                                                                                                                                          |
| `unlock`              | Remove stale Snakemake locks after a cancelled/killed run                                                                                                                                      |
| `propose`             | Run a single `*_propose` rule (optional; not required for the main pipeline)                                                                                                                   |


**Approve aliases:** `qc_review` → `post_qc_review`; `resolution_review` → `s7_clustering`. Parameter keys in `revise` still use internal names (e.g. `s7_clustering.rna.resolution`).

### Local workflow

Planning and QC stages run automatically. Snakemake stops only at the three checkpoints.

```bash
CFG=<run_dir>/deliverables/plan/config/run.yaml

# Confirm execution mode once (required before any compute — run/submit refuse otherwise):
Processing-MuAgent configure-execution --config $CFG --mode local --confirmed-by-user

# Option A — pause at each checkpoint (recommended first time):
Processing-MuAgent run --config $CFG
Processing-MuAgent approve plan_review --config $CFG
Processing-MuAgent run --config $CFG
Processing-MuAgent approve qc_review --config $CFG
Processing-MuAgent run --config $CFG
Processing-MuAgent approve resolution_review --config $CFG
Processing-MuAgent run --config $CFG

# Option B — pre-seed all three checkpoints (unattended Snakemake; you still review outputs):
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
Processing-MuAgent revise s7_clustering s7_clustering.rna.resolution=1.2 --config $CFG
Processing-MuAgent resolution-compare --config $CFG --rna 1.0,1.2 --atac 0.6,0.8
Processing-MuAgent run --config $CFG --no-context
Processing-MuAgent hpc-info
Processing-MuAgent propose post_qc_review --config $CFG
Processing-MuAgent propose s7_clustering --config $CFG
```

## Environment

Recreate the canonical conda env:

```bash
micromamba env create -n grn -f workflow/envs/processing.yaml
micromamba activate grn
pip install -e .
```

