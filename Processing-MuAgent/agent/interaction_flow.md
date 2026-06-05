# Interaction flow — four steps

Canonical entry behaviour for the chat runtime. Each step has three parts:

- **PROMPT_TO_USER** — what the agent says (paraphrase in your own words if natural; preserve the facts).
- **AGENT_ACTIONS** — exact executor CLI invocations to run (no freelancing; these are the only state-changing calls).
- **WHAT_TO_SURFACE_BACK** — what content the agent relays to the user before proceeding.

`$CFG` below is shorthand for the config path passed to every CLI call. After `executor init`, the canonical path is `<run_dir>/deliverables/pre_run/config/run.yaml` — use that value for Steps 2–4.

---

## Step 1 — Declare analysis type

### PROMPT_TO_USER

"Which analysis are you running? I support three options:

1. **scRNA-seq** — single modality, RNA only. Input: 10x h5 / 10x MEX directory / `.h5ad`.
2. **scATAC-seq** — single modality, ATAC only. Input: `fragments.tsv.gz` + matching `.tbi` index.
3. **Paired multiome RNA + ATAC** — one sample, both modalities, cells share barcodes. Inputs: an RNA matrix (formats above) + ATAC fragments. I also handle the *separate* case (RNA + ATAC from independent samples where barcodes don't match) — no cross-modality integration, each modality flows through independently.

For any of these, I do QC → doublet detection → PCA (RNA) + neighbor graph → clustering → UMAP per modality, and then I stop. I don't do integration, annotation, or any downstream analysis.

Before you answer, one other thing I'll need: a run directory — a writable folder where I can put intermediate artifacts, approval checkpoints, and final outputs.

Mandatory approval points depend on the branch:

- All branches: (1) biological context, (2) the full preprocessing plan, (3) **QC review** (after S3, before S6 PCA + neighbor graph), (4) clustering resolution.
- `paired` only: (5) doublet-policy reconciliation at S3 — how to combine the RNA and ATAC detector calls.
- `separate` / `rna_only` / `atac_only`: S3 runs automatically — each modality's doublets are filtered by its own detector with no cross-modal reconciliation; no pause unless you ask for per-stage review.

Which analysis, and where should I write the run?"

### AGENT_ACTIONS

None yet — Step 1 is purely conversational. No executor calls until Step 2.

### WHAT_TO_SURFACE_BACK

After the user answers, confirm back in one line which of the four `workflow_branch` values you inferred (`paired`, `separate`, `rna_only`, `atac_only`) and what the run dir is. Then proceed to Step 2.

See [`stage_prompts/entry.md`](stage_prompts/entry.md) for the canonical Step 1 script and corner cases.

---

## Step 2 — Collect paths + optional biological context

### PROMPT_TO_USER

"Paths I need, based on your declared analysis type:

- **RNA input** (if relevant): full path to `.h5`, MEX directory, or `.h5ad`. Optional **`rna_raw_path`**: companion raw-droplet matrix; when both filtered and raw are supplied, the plan defaults to SoupX for S1a (DecontX when filtered only). **`study_goal`**: `rare_populations` strongly recommends ambient correction; `clustering_inference` defaults to auto but user confirms at plan review. Optional **`s1a_ambient_method`** in run.yaml (`auto`|`none`|`decontx`|`soupx`).
- **ATAC input** (if relevant): full path to `fragments.tsv.gz`. The `.tbi` must sit next to it — I'll fail fast if it doesn't.
- **Genome assembly**: `mm10`, `GRCh38`, etc. Required. I cross-check this against the ATAC fragment chromosome naming, so declare it carefully.
- **Seed** (optional, default 42).

Optional paired-multiome inputs — only collect when relevant:

- **`barcode_translation_path`** (optional): 2-column TSV (`rna_barcode`, `atac_barcode`) mapping cell pairs across whitelists. Required only when GEX and ATAC pipelines used different 10x whitelists (e.g. separate Cell Ranger GEX + Cell Ranger ATAC runs). Without it, S0's diagnostics ladder commits `separate` and the run falls through to the separate branch even if `paired` was declared.
- **`atac_peaks_path`** (optional): BED file of peak intervals. When set, S5 builds the peak-by-cell matrix from these intervals as the highest-priority peak source (ahead of ARC h5 / MACS3 / tile fallback). Spectral embedding is unchanged.
- **`cell_metadata_path`** (optional): per-cell metadata TSV. Must have a `barcode` column for the obs join key. Left-joined into RNA (and ATAC) obs at S8 before the final write. If the file also exposes `rna_barcode`+`atac_barcode` columns, S0's ladder will use those columns as a translation table — so this can stand in for `barcode_translation_path` in a single file.

Biological context is optional but strongly recommended — it shapes QC thresholds and surfaces conflicts early. You can give it to me in any of three forms:

1. **Short chat text** — just tell me: organism, tissue, assay, any related DOIs.
2. **A filled Biological Context Report** — paste the content, or give me a path to an existing markdown file.
3. **A DOI list only** — I'll fetch abstracts and extract what I can.

**Execution environment** — if you haven't said already: should this run **locally** on this machine, or on an **HPC cluster** (PBS Pro or SLURM)? For cluster runs I'll probe available queues/partitions on the login node, suggest a project code or account name where I can detect one, and write `deliverables/pre_run/config/hpc.env` with the `PMA_*` settings you confirm.

Which paths, biological context (if any), and local vs HPC?"

### AGENT_ACTIONS

Once the user answers, in order:

1. Build a minimal `run.yaml` in memory with the fields they supplied:
   ```yaml
   run_dir: <their run dir>
   rna_path: <optional>
   atac_fragments_path: <optional>
   genome_assembly: <required>
   study_goal: clustering_inference     # or rare_populations if user specified
   seed: 42
   ```
   Write it to a temporary path like `<run_dir>/run.yaml.draft` (before init copies it to the canonical location).

2. Invoke `executor init --config <draft-run.yaml>`.
   - This scaffolds `internal/` and `deliverables/{pre_run,checkpoint,post_run}/` and writes the blank biological-context template at `deliverables/pre_run/config/biological_context.md`.
   - From this point forward, use `$CFG = <run_dir>/deliverables/pre_run/config/run.yaml` for every CLI call.

3. Populate biological context (if the user gave any):
   ```python
   from executor import context_mapper
   md = context_mapper.build_report_from_chat(
       organism=..., tissue=..., assay=...,
       dois=[...], notes=...,
   )
   context_mapper.write_report(run_dir, md)
   ```
   If the user gave only a path to an existing filled template, read the file and pass its content to `context_mapper.write_report(run_dir, content)` so it lands at the canonical location. If the user gave nothing, leave the blank template — the Phase 1 gate will block and give them a second chance, or they can explicitly opt out with `executor run --no-context`.

4. **Configure execution mode** (if not already stated in Step 1/2 conversation):
   - If the user said **local** (or gave no preference and you're not on a cluster login node): `executor configure-execution --config $CFG --mode local`.
   - If the user said **HPC** (or you're on a login node with `qsub`/`sbatch` and the dataset is large):
     a. Run `executor hpc-info` and surface the detected scheduler, available queues/partitions, suggested project/account, and any `PMA_*` vars already in the environment.
     b. Ask the user to confirm or override: queue/partition, project/account (if required on their site), optional `PMA_RESOURCES_SCALE`.
     c. Write settings: `executor configure-execution --config $CFG --mode pbs|slurm --pbs-queue ... --pbs-project ...` (or `--slurm-partition` / `--slurm-account`). This records `execution.mode` in `parameters.yaml` and writes `deliverables/pre_run/config/hpc.env`.
   - Do not guess queue/partition/project/account — use `hpc-info` suggestions and ask the user to confirm.

5. Invoke `executor declare-branch <paired|separate|rna_only|atac_only> --config $CFG`.
   - S0 will confirm this assertion against its own pairing detection; if they mismatch, S0 raises with a clear diff — relay it verbatim.

6. **Run planning (P1 → S0 → P2) — S0 execution location depends on configured mode:**
   - **HPC mode (`execution.mode` is `pbs` or `slurm`):** source `hpc.env`, then run S0 on the cluster:
     ```
     source deliverables/pre_run/config/hpc.env
     executor run --config $CFG --executor pbs|slurm --target s0_ingest_execute
     executor run --config $CFG --target p2_plan_execute
     ```
   - **Local mode (`execution.mode` is `local`):** `executor run --config $CFG --target p2_plan_execute` (runs P1 → S0 → P2; ~30s on small inputs).
     - **If S0 fails with a resource error** (OOM, `Killed`, MemoryError, walltime exceeded — use `executor.hpc.looks_like_resource_failure()` on the snakemake stderr):
       1. Configure HPC if not done (`hpc-info` + `configure-execution`; bump `PMA_RESOURCES_SCALE` if appropriate).
       2. `source deliverables/pre_run/config/hpc.env`
       3. `executor run --config $CFG --executor pbs|slurm --target s0_ingest_execute`
       4. `executor run --config $CFG --target p2_plan_execute`
   - **Do not** retry logic errors on the cluster (pairing ambiguous, declared-vs-detected mismatch, missing index). Relay the message; user fixes inputs or branch declaration.
   - After S0 completes, surface `validation_report.json` pairing fields if `paired` was downgraded or ambiguous.

### WHAT_TO_SURFACE_BACK

- Confirm `executor init` wrote `deliverables/pre_run/config/run.yaml` and `biological_context.md`.
- If context was supplied, confirm it was written (`deliverables/pre_run/config/biological_context.md`).
- If HPC mode was configured, confirm `execution.mode` and the path to `deliverables/pre_run/config/hpc.env`; remind the user to `source` it before cluster submit/resume.
- After `executor run --target p2_plan_execute`, surface `deliverables/pre_run/summary/context_summary.md` if populated (conflicts or inferred values). Do not paraphrase — paste the markdown back.

See [`stage_prompts/inputs_intake.md`](stage_prompts/inputs_intake.md) for the canonical Step 2 script and the per-context-form handling details.

---

## Step 3 — Confirm the plan

### PROMPT_TO_USER

"Here's the preprocessing plan. Nothing heavy runs until you approve it — this is the single point where cross-stage inconsistencies are cheapest to catch. Review and tell me: **approve as-is**, **revise one or more parameters**, or **abort**.

[... plan_review.md content relayed verbatim ...]"

### AGENT_ACTIONS

1. Invoke `executor plan-review --config $CFG`.
   - This re-renders (and writes) `deliverables/pre_run/summary/plan_review.md` — summary section (8 decision items) plus a full parameter appendix.
   - The same content also lives at that path if the user wants to open it directly.

2. On user decision:
   - **Approve** → `executor approve plan_review --config $CFG --note "approved after review"`.
   - **Revise a parameter** → `executor revise <stage> <key>=<value> --config $CFG --rationale "<user's reason>"`. Stage is re-set to awaiting_approval; ask if more revisions are needed before re-approving.
   - **Abort** → stop. Tell the user the run dir is intact; they can resume later by re-invoking you on the same config.

### WHAT_TO_SURFACE_BACK

The **Summary** section of `plan_review.md`, verbatim (the appendix is reference detail — surface decision points and any `?` flags from the summary). Do not paraphrase values.

---

## Step 4 — Run with checkpoints

### PROMPT_TO_USER

(At each mandatory pause, surface the stage's proposal content and ask for approval or revision. Pauses are branch-aware:)

- **`p1_context`** (all branches): biological context extraction + conflict resolution. Already handled in Step 2 flow in most cases, but if the user skipped context in Step 2, P1 will stop here.
- **`plan_review`** (all branches): covered in Step 3 — checkpoint **#1**.
- **`post_qc_review`** (all branches): QC review checkpoint **#2** between S3 and S4/S5. Generates QC figures and `checkpoint/qc_review/qc_review.md` (S1–S3 metrics; on **paired**, includes S3 union doublet policy for confirmation). Point the user at `deliverables/checkpoint/qc_review/`. They may revise S1/S2 thresholds and re-run affected stages before approving. On `separate` / single-modality branches, no cross-modal doublet policy applies.
- **`s7_clustering`** (all branches): resolution review checkpoint **#3**. Review `checkpoint/resolution_review/`. **Separate / single-modality:** resolutions set **final** labels in processed outputs. **Paired:** **diagnostic** per-modality labels for UMAP only.
- **`s3_doublets`**: not a separate user checkpoint — runs before QC review; policy is confirmed at checkpoint **#2** on paired runs. Auto-approve unless the user asked for stage-by-stage review.

(For other stages, auto-approve silently unless the user asked for per-stage review.)

"I'm at the **[stage]** checkpoint. Here's what the pipeline proposes:

[... proposal yaml + summary markdown content, verbatim ...]

Approve, revise, or abort?"

### AGENT_ACTIONS

1. Invoke `executor run --config $CFG` (no `--auto-approve`) for **local** mode, or follow the HPC table below when `execution.mode` is `pbs`/`slurm` (read from `parameters.yaml` or ask if missing).
   - Snakemake runs every stage whose approval sentinel exists, stopping at the first stage whose `.approved` is missing.
   - Loop:
     a. Run `executor status --config $CFG` to see which stage is currently `awaiting_approval`.
     b. Read `<run_dir>/internal/proposals/<stage>.yaml` — the structured proposal.
     c. If the stage has a linked summary in `deliverables/checkpoint/resolution_review/` (e.g., `resolution_summary.md` for s7), read that too and surface both.
     d. Based on user decision:
        - Approve → `executor approve <stage> --config $CFG`.
        - Revise → see **QC threshold revision** below if the current stage is `post_qc_review`; otherwise `executor revise <stage> <key>=<value> --config $CFG`; re-surface the updated proposal; loop.
     e. Re-invoke the appropriate run command (see HPC section below). Continue until `manifest` completes.

**QC threshold revision at `post_qc_review` (detailed procedure):**

When the user asks to adjust an S1, S2, or S3 parameter at the QC review checkpoint, execute these steps in order — do not skip any:

1. **Update parameters.** For each changed parameter run:
   ```
   executor revise <stage> <stage>.<param>=<value> --config $CFG --rationale "<user's reason>"
   ```
   This writes the new value to `parameters.yaml` and marks the stage `awaiting_approval`.

2. **Delete stale artifacts** so Snakemake re-runs the affected stages:
   - **S1 revised:** `internal/artifacts/s1_rna_qc/rna_qc.h5ad` — plus all S3 artifacts listed below.
   - **S2 revised:** `internal/artifacts/s2_atac_qc/atac_qc.h5ad`, `atac_snap.h5ad`, `qc_summary.json`. Keep `atac_fragments_cbf_chrnorm.tsv.gz*` (expensive to regenerate) — plus all S3 artifacts listed below.
   - **S3 revised:** `internal/artifacts/s3_doublets/rna_post_doublet.h5ad`, `atac_post_doublet.h5ad`, `calls.parquet`, `joint_barcodes.txt`, `overlap_summary.json`.
   - Any revision at S1 or S2 also invalidates S3 — always delete S3 artifacts too.

3. **Approve revised stages.** For each stage whose artifacts were deleted:
   ```
   executor approve <stage> --config $CFG
   ```
   Do **not** pass `--auto-approve` to `submit` (it refreshes timestamps on all sentinels and can trigger spurious re-runs of already-complete stages).

4. **Submit.** Resubmit the pipeline (HPC mode):
   ```
   source deliverables/pre_run/config/hpc.env
   executor submit --config $CFG --executor pbs|slurm
   ```
   Monitor with `executor hpc-status --watch --config $CFG`.

5. **Regenerate QC reports.** The submit target is `s3_doublets_execute`; the head job exits after S3 without running the local propose rule. After S3 completes, run the propose rule explicitly on the login node:
   ```
   executor propose post_qc_review --config $CFG
   ```

6. **Surface the updated report.** Read and relay `deliverables/checkpoint/qc_review/qc_review_<run_name>.md` verbatim. Ask the user to approve or revise again.

2. When `manifest` finishes:
   - Read `deliverables/post_run/run_manifest.json` and extract `workflow_branch`, `outputs`.
   - Point the user at `deliverables/post_run/qc_summary.md`, the UMAP figures in `deliverables/post_run/`, and the handoff artifact `run_manifest.json`.

### WHAT_TO_SURFACE_BACK

- At each pause: the full proposal yaml content + any linked summary markdown.
- At completion: the manifest's `outputs` keys + a one-line sign-off ("Run complete. Outputs at `deliverables/post_run/`. I stop here — integration/annotation is out of scope.").

---

## Running on HPC (PBS Pro or SLURM)

The four-step flow above is unchanged on a cluster. Only the underlying execution
model differs: heavy `_execute` rules dispatch to scheduler jobs, while every
`*_propose` rule and the light planning stages (`p1_context`, `p2_plan`,
`plan_review`) plus `manifest` are declared `localrules` and run on the orchestrator
host (login node in interactive mode; the head-job in headless mode).

**S0 ingest** runs on the cluster directly when `execution.mode` is `pbs` or
`slurm`; runs locally when `execution.mode` is `local` (with cluster retry on
resource failure). Pairing-detection conflicts are always resolved interactively
on the login node after S0 finishes.

### Intake (Step 2) — ask before P1 runs

If the user has not said whether to run locally or on HPC, ask alongside the
biological-context question. When they choose HPC:

1. Run `executor hpc-info` on the login node.
2. Surface detected scheduler, available queues/partitions, suggested project or
   account (from env vars or recent jobs), and any `PMA_*` vars already set.
3. Ask the user to confirm or override queue/partition, project/account, and
   optional `PMA_RESOURCES_SCALE`.
4. Run `executor configure-execution --config $CFG --mode pbs|slurm ...` to write:
   - `deliverables/pre_run/config/hpc.env` — shell snippet sourced by runner scripts
   - `deliverables/pre_run/config/site.config` — YAML platform description consumed by Execution-MuAgent
   - Records `execution.mode` in `parameters.yaml`.

Do not invent queue or account names — probe with `hpc-info`, suggest, and wait
for confirmation.

### HPC run phases (after plan review)

| Step | Stages | Executes on | You |
|------|--------|-------------|-----|
| Planning | P1 → P2 | Login node | — |
| S0 ingest | S0 | Cluster (HPC mode) / Login node (local mode) | — |
| Checkpoint **#1** | plan_review | Login node | Review plan |
| QC | S1a → S3 | Cluster | — |
| Checkpoint **#2** | post_qc_review | — | Review QC |
| PCA + neighbors + clustering | S4 → S7 (sweep) | Cluster | — |
| Checkpoint **#3** | s7_clustering | — | Review resolution |
| Finish | S7 (labels) → S8 → manifest | Cluster | — |

After plan review approval, `source deliverables/pre_run/config/hpc.env`, then:

- **QC batch:** `executor submit --config $CFG --executor pbs|slurm --auto-approve --auto-approve-except post_qc_review --auto-approve-except s7_clustering`
- **After QC approval:** same submit with `--auto-approve-except s7_clustering` only
- **After resolution approval:** `executor submit --config $CFG --executor pbs|slurm`

`executor submit` is a hard dependency on Execution-MuAgent — it writes `internal/stage_meta/head_job.yaml`, then starts `Execution-MuAgent execute-spec` as a **background supervision daemon**. The daemon submits the head-job, records the job ID to `execution_manifest.jsonl`, and then runs the watch loop (stall detection + kill-on-hang) for the full lifetime of the job. `submit` returns within ~90 seconds, as soon as the job ID is confirmed. If Execution-MuAgent is absent, `submit` fails loudly.

The daemon survives SSH disconnects on most clusters. On sites with `KillUserProcesses=yes`, remind the user to run inside `tmux` or `screen`. The daemon writes its output to `internal/hpc_monitor/monitor_<timestamp>.log` (with a `monitor.log` symlink to the latest) and removes `monitor.pid` when it exits.

Poll job health with `executor hpc-status --watch --config $CFG`. Findings and hang reports appear in `internal/hpc_monitor/latest_report.md`.

### How Claude surfaces things during the HPC flow

- **At Step 2**, if execution mode is unknown, ask local vs HPC before `executor run --target p2_plan_execute`. Run `hpc-info` and walk through `hpc.env` setup when HPC is chosen.
- **At Step 1–4 otherwise**, behaviour matches local mode until plan review is approved.
- **When the user runs `submit`**, surface the printed PBS/SLURM job ID, the supervision daemon PID, and the log path. Remind them that the daemon is running in the background and they can follow job health with `Processing-MuAgent hpc-status --watch --config $CFG`.
- **When `status --watch` shows `s7_clustering awaiting_approval`**, surface the path to `resolution_review.html` (the
  primary review artifact) AND `resolution_review.ipynb` (for power users who
  want to re-cluster at custom resolutions interactively). Paste the contents
  of `resolution_summary.md` verbatim, plus the proposal's `review_artifacts`
  block.
- **On approve/revise**, run the appropriate CLI as today, then `submit` again
  (or `run --executor pbs|slurm` if the user prefers foreground on the login node).

### Site variables (user confirms after `hpc-info`)

```bash
# PBS Pro example:
export PMA_PBS_QUEUE=<your_queue_name>
export PMA_PBS_PROJECT=<your_project_code>

# SLURM example:
export PMA_SLURM_PARTITION=<your_partition_name>
export PMA_SLURM_ACCOUNT=<your_account_name>

# Optional — scale per-rule memory and walltime (default is 1):
export PMA_RESOURCES_SCALE=2
```

These are written to `deliverables/pre_run/config/hpc.env` by `configure-execution`.

---

## When things go wrong

- **S0 auto-downgrades `paired → separate`** (not an error — this is the new diagnostics-ladder behaviour). User declared `paired` but none of the ladder rungs (direct overlap, suffix-normalized, `barcode_translation_path`, `cell_metadata_path` with `atac_barcode` column) validated cell-level pairing. The run continues on the `separate` branch. Surface `internal/artifacts/s0_ingest/validation_report.json` verbatim — its `pairing.downgrade_reason` field has the specific reason. Ask the user whether to (a) proceed on `separate`, (b) supply a `barcode_translation_path` and rerun S0, or (c) abort and fix inputs upstream (e.g. rerun `cellranger-arc count` for a combined ARC matrix).
- **S0 raises "declared=... conflicts with detected=..."** — the user's `executor declare-branch` doesn't match what S0 detected, AND the declaration is single-modality (`rna_only`/`atac_only`). Single-modality conflicts still raise hard because they signal a data-hygiene problem. Relay the raised message. Ask the user to either correct the declaration or correct the config (drop the unwanted modality).
- **S0 raises "pairing is ambiguous"** — RNA+ATAC Jaccard overlap is between 30% and 80% after normalization/subset checks, and the user did not declare a branch (or declared one that doesn't resolve the ambiguity). Ask the user: are these paired or separate? Based on answer, run `executor declare-branch <paired|separate>` and re-run; for `paired`, supply `barcode_translation_path` if barcode whitelists differ. Don't auto-pick.
- **S3 raises "paired-branch joint barcode intersection is empty"** — S0 committed `paired` but no cell survived both modalities' QC + doublet removal. This usually means QC thresholds were too aggressive. Surface the message; ask the user to revise S1/S2 thresholds via `executor revise s1_rna_qc ...` or `executor revise s2_atac_qc ...`. If the pairing decision used `pairing.translation_table`, also check that the translation table actually covers the QC-surviving cell set.
- **Phase 1 gate raises "biological_context.md is empty"** — the user didn't give context and didn't opt out. Ask for context OR offer `executor run --config $CFG --no-context` as the explicit opt-out.
- **S0 fails with OOM / Killed / walltime on the login node** — not a logic error. Configure HPC if needed, `source hpc.env`, retry `executor run --config $CFG --executor pbs|slurm --target s0_ingest_execute`, then resume `executor run --config $CFG --target p2_plan_execute`. Consider raising `PMA_RESOURCES_SCALE`. Tell the user what you did.
- **A stage execute fails at runtime** — relay the traceback from snakemake. Do not retry silently; root-cause first. If the user insists on retry, use `executor run --config $CFG --target <stage>_execute`.
- **Execution-MuAgent reports `submit_rejected_policy`** — the scheduler rejected the job as a policy error (invalid partition, account, or walltime over the site limit). Read `internal/hpc_monitor/latest_report.md` for the scheduler's exact message. Tell the user which field to correct: partition/account via `executor configure-execution --mode <scheduler> ...` (rewrites `site.config`), or walltime by reducing `PMA_RESOURCES_SCALE`. Then `executor submit` again.
- **Per-stage specs not written** — specs are written automatically by `executor plan-review`. If `internal/specs/` is missing or empty, re-run `executor plan-review --config $CFG`. Specs are internal state; do not surface them to the user unless asked.
- **`hpc-status` shows "Supervisor: not running" alongside a RUNNING or PENDING scheduler state** — the supervision daemon has died but the cluster job is still active. Without the daemon, stalled jobs will not be auto-cancelled. Restart the daemon: `executor supervisor-restart --config $CFG`. This resumes the full watch loop (stall detection, kill-on-hang) against the already-running job without resubmitting. Tell the user what happened and what you did.
- **Supervision daemon crashes on a site with KillUserProcesses=yes** — when the user's SSH session ends, systemd kills all their processes including the daemon. The cluster job keeps running, but protection is gone. For the current run, tell them to use `supervisor-restart` as soon as they reconnect. Going forward, suggest running `submit` inside a `tmux` or `screen` session on that cluster.
