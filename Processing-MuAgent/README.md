# Processing-MuAgent

Multiome (scRNA-seq + scATAC-seq) preprocessing subagent. Takes filtered or raw 10x Genomics multiome outputs and performs QC, PCA (RNA) + neighbor graphs, clustering, and UMAP **per modality**, then **stops** before integration.

Supported workflow branches: `paired`, `separate`, `rna_only`, `atac_only`. Declare the branch up front with `Processing-MuAgent declare-branch`.

## Pipeline overview

**Stage order** (Snakemake DAG — ingest validation must finish before preprocessing plan assembly, because P2 reads `internal/artifacts/s0_ingest/validation_report.json`):

```
  P1 context extraction → S0 ingest validation → P2 preprocessing plan → (CHECKPOINT 1) plan_review
  → S1a ambient RNA correction → S1 RNA QC → S2 ATAC QC → S3 doublets → (CHECKPOINT 2) post_qc_review
  → S4 RNA normalization + HVG → S5 ATAC spectral embedding → S6 PCA (RNA) + neighbor graph
  → S7 clustering + (CHECKPOINT 3) resolution_review → S8 UMAP → outputs
```

Automated processing stages (`p1_context`, `s0_ingest`, `p2_plan`, `s1a`–`S3`, `S4`–`S6`, `s8_umap`) run from **artifact dependencies** and the three checkpoint boundaries below — they do **not** require per-stage `approve` calls. Optional `<stage>_propose` rules still exist for inspection or debugging, but they are not on the main execution path.

### User checkpoints (3)

Three deliberate pauses where you review deliverables and decide before heavy downstream work continues. Everything else runs automatically once upstream artifacts exist and the relevant checkpoint is approved.

| # | CLI name | Internal stage | When | What you decide |
|---|----------|----------------|------|-----------------|
| **1** | **Plan review** | `plan_review` | After S0 + P2, before S1 | Approve the preprocessing plan (`pre_run/summary/plan_review.md`) |
| **2** | **QC review** | `post_qc_review` | After S3, before S4/S5 | Inspect QC figures + `checkpoint/qc_review/qc_review.md`; revise **RNA/ATAC quality-filter thresholds** and re-run if needed; on **paired** multiome, confirm the **union doublet removal policy** |
| **3** | **Resolution review** | `s7_clustering` | After S6, before S8 | Choose Leiden resolution per modality from sweep metrics (`checkpoint/resolution_review/`). **Separate / single-modality:** sets **final** cluster labels. **Paired:** **diagnostic** per-modality labels for UMAP only (not joint embedding) |

## Workflow stages

### Planning (pre-QC)

- **P1 Context extraction** — Biological Context Report (organism, tissue, assay, DOIs) plus DOI-based prior-analysis extraction.
- **S0 Ingest** — Accepts Cell Ranger **filtered** and **raw** matrices, auto-detecting RNA and ATAC formats (see tables below) and validating fragments files. Performs a **diagnostic barcode check for paired multiome**: S0 checks for direct barcode matches, then tries matching after removing suffixes. If those don't match, it looks for a `barcode_translation_path` or `cell_metadata_path` provided by the user. No barcode intersection is performed here. If pairing can’t be confirmed, S0 downgrades the workflow from `paired` to `separate` and logs the reason in `validation_report.json`. S0 runs after P1 and before P2; its validation report feeds into the preprocessing plan. When HPC is configured (`execution.mode = pbs | slurm`), S0 runs as a **cluster job** directly. In local mode, S0 runs on the login node and is retried on the cluster only if it hits a resource limit (OOM or walltime).

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

- **S1a Ambient RNA correction** — Default `method=auto` on RNA branches (SoupX if raw+filtered exist, else DecontX). Omitted on `atac_only`. Whether to run is confirmed by user at **plan review (#1)** depending on the study goal, inputs, and sample context (see [10x ambient RNA guide](https://www.10xgenomics.com/analysis-guides/introduction-to-ambient-rna-correction)).
- **S1 RNA QC** — MAD-derived thresholds on `total_counts` / `n_genes` / `pct_counts_mt` plus a `pct_counts_ribo` ceiling, computed on decontaminated counts from S1a. Writes pre/post QC violin figures to `deliverables/checkpoint/qc_review/`.
- **S2 ATAC QC** — Four per-cell filters applied in sequence:
  1. **n_fragments** — MAD-based bounds on log-scale fragment count (with an absolute floor).
  2. **TSS enrichment** — min/max bounds on SnapATAC2's TSS score.
  3. **Nucleosome signal** — upper bound on Signac-style `mono / nucleosome_free` ratio.
  4. **FRiP** — Fraction of Reads in Peaks (`frip_min`; default 0.15). S2 acquires a peak set using the same priority order as S5 feature export (see below), calls the peak × cell matrix via SnapATAC2, and filters cells with FRiP below threshold. The peak BED written here is reused by S5 so no redundant peak calling occurs. FRiP filtering is silently skipped when no peak source is available.

  Writes fragment-size distribution figures to `deliverables/checkpoint/qc_review/`.
- **S3 Doublets** — Per-modality doublet detection, then branch-specific reconciliation:
  - **RNA:** Scrublet (sparse-CSR input; `expected_doublet_rate ≈ 0.0008 × n_cells`, capped at 10%).
  - **ATAC:** SnapATAC2 scrublet (thresholds configurable in the preprocessing plan).
  - **separate / single-modality branches:** Each modality is filtered independently by its own detector; per-modality calls are saved in `calls.parquet`.
  - **paired branch:** Also performs joint barcode alignment after doublet removal; the union doublet policy is confirmed at the **QC review checkpoint** (`checkpoint/qc_review/qc_review.md`).
- **post_qc_review** — **QC review checkpoint (#2).** Propose-only gate between S3 and S6 PCA (RNA) + neighbor graph. Generates doublet histograms, a cell-count waterfall (with counts labelled on bars), and `checkpoint/qc_review/qc_review.md` — a plain-language summary of what each filter step did (MAD outlier bounds, MT/ribo ceilings, TSS enrichment, nucleosome signal, FRiP, union doublet policy). Revise quality-filter thresholds and re-run affected stages before approving.
- **S4 RNA norm + HVG** — Log-normalize (`target_sum=1e4`) + HVG selection (`seurat_v3` on counts).
- **S5 ATAC spectral embedding and peak matrix export** — SnapATAC2 tile matrix (`bin_size=500`, unified with S3) → feature selection → `snap.tl.spectral` (Laplacian eigenmaps with IDF feature weights; not classical TF-IDF + SVD LSI). In parallel, exports a feature (cell-by-feature) matrix using this priority order for the peak coordinates:
  0. **User-supplied peaks** — `atac_peaks_path` in `run.yaml` → SnapATAC2 `make_peak_matrix` (`user_peaks` mode).
  1. **ARC peak matrix** — pre-called peaks from a combined Cell Ranger ARC `.h5` detected at S0 (`arc_h5` mode).
  2. **S2 pre-called peaks** — BED file written by S2 ATAC QC (MACS3 or ARC-derived) reused here; no redundant peak calling (`s2_peaks_macs3` / `s2_peaks_arc` mode).
  3. **Tile-matrix fallback** — verified SnapATAC2 tile matrix (`tile_matrix_fallback` mode), used only when no peak source is available.

  Spectral embedding in `obsm['X_spectral']` (with `X_lsi` as a backward-compat alias) is always computed from the tile matrix regardless of peak-export mode. When `drop_first=True`, the first component is removed before S6–S8.
- **S6 PCA (RNA) + neighbor graph** (`s6_neighbors`) — **RNA:** optional `sc.pp.scale`, then PCA; `n_pcs` from a chord-distance elbow on explained variance, capped at `rna_n_pcs_max`; nearest-neighbors on PCA space. **ATAC:** KNN graph on the S5 spectral embedding (`X_spectral` via `snap.pp.knn`). Artifact: `internal/artifacts/s6_neighbors/rna_neighbors.h5ad`.
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

# Write site.config (source of truth) and derived hpc.env:
# PBS Pro:
Processing-MuAgent configure-execution --config $CFG --mode pbs \
  --pbs-queue <queue> --pbs-project <project>

# SLURM:
Processing-MuAgent configure-execution --config $CFG --mode slurm \
  --slurm-partition <partition> --slurm-account <account>
```

This writes:
- `deliverables/pre_run/config/site.config` — YAML platform description (consumed by Execution-MuAgent)
- `deliverables/pre_run/config/hpc.env` — shell snippet generated from site.config; source before `submit`

For larger datasets increase `--resources-scale` (e.g. `2` for ~30k cells, `4` for ~100k). Per-stage CPU, memory, and walltime defaults live in `workflow/resources.smk`; OOM-killed jobs are retried once at double memory.

### Requirements

- **SLURM:** Uses Snakemake's generic cluster executor plus `sbatch`/`sacct`/`squeue`; every heavy rule is submitted as a separate `sbatch` job, even when the Snakemake head-job itself runs under SLURM. Child jobscripts are sanitized before submission so Snakemake does not re-enable `storage-local-copies` on shared NFS (a common cause of jobs finishing Python quickly then hanging at `Storing output in storage.`).
- **Clean stage exit:** after writing its outputs, every cluster `<stage>_execute` rule calls `executor.cluster_exit.finalize_cluster_exit()` (`gc.collect()` + `os._exit(0)` when `SLURM_JOB_ID`/`PBS_JOBID` is set). This terminates lingering h5py/HDF5 background threads that would otherwise keep the child process alive after the output is complete — the previous cause of jobs stuck RUNNING with `slurmstepd: error: Pid still in cpuset cgroup`. It is a no-op in local mode (job-id env vars unset).

### How the HPC run proceeds

| Step | Stages | Executes on | You |
|------|--------|-------------|-----|
| Planning | P1 → P2 | Login node (default), or `srun` on a compute node if the login node memory is limited| — |
| S0 ingest | S0 | Cluster (HPC mode) / Login node (local mode) | — |
| Checkpoint **#1** | plan_review | Login node | Review plan |
| QC | S1a → S1 → S2 → S3 | Cluster | — |
| Checkpoint **#2** | post_qc_review | — | Review QC |
| PCA + neighbors + clustering | S4 → S5 → S6 → S7 (sweep) | Cluster | — |
| Checkpoint **#3** | s7_clustering | — | Review resolution |
| Finish | S7 (labels) → S8 → manifest | Cluster | — |

**S0 execution mode:** in HPC mode (`execution.mode = pbs | slurm`), the agent runs S0 on the cluster directly (after sourcing `hpc.env`). In local mode, S0 runs on the login node; on OOM or walltime failure the agent configures HPC (if needed) and retries `s0_ingest_execute` on the cluster before continuing with P2.

Each heavy `_execute` stage runs as its own scheduler job. Only the three checkpoints above require `approve`. Findings and hang reports are written to `internal/hpc_monitor/latest_report.md` by Execution-MuAgent.

### Submit workflow

After checkpoint **#1** (`plan_review`), source `deliverables/pre_run/config/hpc.env`, then use `Processing-MuAgent submit` (not `run`) to dispatch the Snakemake head-job. **`submit` auto-infers the Snakemake target** from checkpoint state — you do not need to pick `post_qc_review_propose`, `s7_clustering_propose`, or `all` manually. After each approval, run `submit` again and it stops at the next gate:

| Checkpoint state | Inferred target | Runs through |
|------------------|-----------------|--------------|
| `post_qc_review` not approved | `post_qc_review_propose` | S1a → S3 + QC summary, then pauses |
| `s7_clustering` not approved | `s7_clustering_propose` | S4 → S6 + resolution sweep, then pauses |
| Both approved | `all` | S7 labels → S8 → manifest |

Override with `--target <name>` only when debugging.

If a previous run was cancelled or killed, stale Snakemake locks under `<run_dir>/internal/snakemake/` can make the next submit fail immediately with `LockException`. Recover with:

```bash
Processing-MuAgent unlock --config $CFG
# or on submit:
Processing-MuAgent submit --config $CFG --executor slurm --unlock-stale-locks
```

```bash
CFG=<run_dir>/deliverables/pre_run/config/run.yaml

# Configure HPC settings (writes site.config + derived hpc.env):
Processing-MuAgent configure-execution --config $CFG --mode slurm \
  --slurm-partition cpu-medium --slurm-account mylab

source <run_dir>/deliverables/pre_run/config/hpc.env

# Run planning + plan review on login node, then submit heavy batch:
Processing-MuAgent run --config $CFG --target p2_plan_execute
Processing-MuAgent plan-review --config $CFG  # also writes internal/stage_meta/
Processing-MuAgent approve plan_review --config $CFG

# First heavy batch (stops at QC review):
Processing-MuAgent submit --config $CFG --executor slurm

# After QC review:
Processing-MuAgent approve qc_review --config $CFG
Processing-MuAgent submit --config $CFG --executor slurm

# After resolution review:
Processing-MuAgent approve resolution_review --config $CFG
Processing-MuAgent submit --config $CFG --executor slurm
```

Each `submit` call starts a background supervision daemon that watches for stalls and cancels hung jobs — it keeps running after `submit` returns. Poll job health and per-step state with `Processing-MuAgent hpc-status --watch --config $CFG`. For an unattended batch, pre-seed all three checkpoints with `--auto-approve` on `run` or `submit`. To keep specific gates interactive, add `--auto-approve-except qc_review` (repeatable; accepts aliases or internal names).

### Execution-MuAgent integration

`Processing-MuAgent submit` delegates the full submission lifecycle to the sibling `Execution-MuAgent` package (a hard dependency — the command fails loudly if it is absent):

1. `plan-review` writes per-stage metadata YAMLs to `internal/stage_meta/` (resources, I/O paths, `progress_timeout_hint`, science description). `progress_timeout_hint` values come from `resources.smk` — that is the single source of truth; `--stale-minutes 90` is the fallback default in Execution-MuAgent when no hint is present.
2. `configure-execution` writes `site.config` (the single platform source of truth). `hpc.env` is generated from it — the two cannot drift.
3. `submit` writes `internal/stage_meta/head_job.yaml`, then starts `Execution-MuAgent execute-spec` as a **background supervision daemon**. The daemon submits the head-job, writes `execution_manifest.jsonl`, and then watches the job for the rest of its lifetime. `submit` returns within ~90 seconds (as soon as the job ID appears in the manifest). Snakemake submits per-stage child jobs from within the head-job.

**The supervision daemon is a safety layer, not optional monitoring.** It is the process that detects stalls, gathers evidence, and cancels hung jobs. If it dies, the cluster job continues to run unprotected — stalled jobs will not be auto-cancelled. `hpc-status` shows whether the daemon is alive and warns you if it isn't.

**Daemon lifecycle.** The daemon runs as a detached background process (`start_new_session=True`) that survives SSH disconnects on most systems. Exception: sites configured with `KillUserProcesses=yes` (a systemd setting) will kill all your processes on logout regardless. On those systems, wrap your session in `tmux` or `screen` before running `submit`. The daemon's output goes to `internal/hpc_monitor/monitor_<timestamp>.log` (with `monitor.log` as a symlink to the latest). The daemon removes `monitor.pid` when it exits normally.

**If the daemon dies mid-run** (SSH drop on a KillUserProcesses site, OOM, unexpected crash) while the cluster job is still active, restart it without resubmitting:

```bash
Processing-MuAgent supervisor-restart --config $CFG
```

This picks up from `latest_submission.json` and resumes the full watch loop (stall detection + kill-on-hang protection) against the already-running job.

**Job monitoring.** The daemon checks in every 4.5 minutes (270 s). It counts how many checks have passed with no sign of progress — no file updates, no log growth. After enough quiet checks (determined by each stage's declared `progress_timeout_hint`), it switches to investigation mode: it looks at CPU and memory usage, tests filesystem responsiveness, and checks for error keywords in logs. A stall signal is suspicion, never a verdict — the daemon only cancels a job when the evidence is conclusive. It always kills children first, then the head job, so Snakemake can resubmit cleanly.

**Per-step output verification.** On every check the daemon verifies each stage's declared outputs as they appear — it actually opens the file (h5ad, parquet, JSON), not just checks for non-empty size. A `stage_output_verified` finding is written when a stage completes cleanly. When the head-job reaches COMPLETED, the same check runs over every stage spec; any missing or corrupt output is reported as `output_missing`.

**When a job is cancelled or fails.** The daemon writes a plain-English diagnosis to `internal/hpc_monitor/latest_report.md` and a machine-readable snapshot to `latest_snapshot.json`, then stops. It never contacts a human and never resubmits. Read the report with `Execution-MuAgent report --run-dir <run_dir>`, fix the underlying problem, and re-run `submit`.

If Execution-MuAgent reports a **policy rejection** (`submit_rejected_policy` in `latest_report.md`), the scheduler refused the job because of an invalid partition, account, or walltime. Re-run `configure-execution` with corrected settings and resubmit.

All cluster submission and monitoring goes through `Processing-MuAgent submit` → `Execution-MuAgent execute-spec`. There is no manual-submission path — Execution-MuAgent runs only via Processing-MuAgent.


## Run directory layout

Per-run state lives under `run_dir` from your config — never inside the source tree.

```
<run_dir>/
  deliverables/
    pre_run/
      config/
        run.yaml                  ← canonical config (use this for all CLI calls)
        biological_context.md     ← Biological Context Report
        hpc.env                   ← source before submit/run on cluster
        site.config               ← YAML platform description (consumed by Execution-MuAgent)
      summary/
        context_summary.md        ← P1 output
        plan_review.md            ← plan review gate (summary + parameter appendix)
    checkpoint/
      qc_review/                  ← QC review checkpoint (#2): figures + qc_review.md
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

Commands below use `$CFG` = `<run_dir>/deliverables/pre_run/config/run.yaml` (written by `init`).

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
CFG=<run_dir>/deliverables/pre_run/config/run.yaml
Processing-MuAgent declare-branch paired --config $CFG   # paired | separate | rna_only | atac_only
```

`init` creates `<run_dir>/`, copies config to `deliverables/pre_run/config/run.yaml`, and writes the Biological Context Report template at `deliverables/pre_run/config/biological_context.md`.

### Command reference

| Command | Purpose |
|---------|---------|
| `init` | Create run directory scaffold |
| `declare-branch` | Record workflow branch in `parameters.yaml` |
| `configure-execution` | Set `execution.mode`; write `site.config` (platform source of truth) and derived `hpc.env` |
| `hpc-info` | Probe PBS/SLURM queues, partitions, accounts on the login node |
| `run` | Foreground Snakemake (`--executor local\|pbs\|slurm`) |
| `submit` | Submit head-job via Execution-MuAgent (hard dependency); starts background supervision daemon; infers phase target |
| `supervisor-restart` | Restart the supervision daemon without resubmitting — use when the daemon died mid-run (SSH drop, site reboot, OOM) but the cluster job is still active |
| `status` | Per-step pipeline state (S1a–S8 + review gates); `--watch` polls until actionable |
| `hpc-status` | Job health (HEALTHY/SUSPECT/…), supervisor liveness, and per-step state; `--watch` polls; warns if supervision is offline while the cluster job is still running |
| `approve` | Write `internal/checkpoints/<stage>.approved` (human checkpoints only) |
| `plan-review` | Render `plan_review.md`; also writes per-stage metadata to `internal/stage_meta/` |
| `revise` | Update one parameter in `parameters.yaml` and reset a checkpoint to awaiting |
| `resolution-compare` | Side-by-side resolution comparison figures (optional) |
| `unlock` | Remove stale Snakemake locks after a cancelled/killed run |
| `propose` | Run a single `*_propose` rule (optional; not required for the main pipeline) |

**Approve aliases:** `qc_review` → `post_qc_review`; `resolution_review` → `s7_clustering`. Parameter keys in `revise` still use internal names (e.g. `s7_clustering.rna.resolution`).

### Local workflow

Planning and QC stages run automatically. Snakemake stops only at the three checkpoints.

```bash
CFG=<run_dir>/deliverables/pre_run/config/run.yaml

# Option A — pause at each checkpoint (recommended first time):
Processing-MuAgent run --config $CFG --executor local
Processing-MuAgent approve plan_review --config $CFG
Processing-MuAgent run --config $CFG --executor local
Processing-MuAgent approve qc_review --config $CFG
Processing-MuAgent run --config $CFG --executor local
Processing-MuAgent approve resolution_review --config $CFG
Processing-MuAgent run --config $CFG --executor local

# Option B — pre-seed all three checkpoints (unattended Snakemake; you still review outputs):
Processing-MuAgent run --config $CFG --executor local --auto-approve

# Option C — unattended except one gate (example: keep QC review interactive):
Processing-MuAgent run --config $CFG --executor local --auto-approve --auto-approve-except qc_review

Processing-MuAgent status --watch --config $CFG
```

`run` requires a filled Biological Context Report unless you pass `--no-context`. After `revise`, approve the affected checkpoint again before resuming.

### HPC workflow

After checkpoint **#1**, use `submit` instead of foreground `run`. See **Running on HPC → Submit workflow** for the resume loop, `unlock`, and `--unlock-stale-locks`.

```bash
source <run_dir>/deliverables/pre_run/config/hpc.env
Processing-MuAgent submit --config $CFG --executor slurm
# submit returns within ~90 s; the supervision daemon keeps running in the background.
Processing-MuAgent hpc-status --watch --config $CFG   # shows job health, supervisor liveness, per-step state

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