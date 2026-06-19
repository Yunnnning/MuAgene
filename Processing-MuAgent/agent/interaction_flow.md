# Interaction flow — four steps

Canonical entry behaviour for the chat runtime. Each step has three parts:

- **PROMPT_TO_USER** — what the agent says (paraphrase in your own words if natural; preserve the facts).
- **AGENT_ACTIONS** — exact executor CLI invocations to run (no freelancing; these are the only state-changing calls).
- **WHAT_TO_SURFACE_BACK** — what content the agent relays to the user before proceeding.

`$CFG` below is shorthand for the config path passed to every CLI call. After `executor init`, the canonical path is `<run_dir>/deliverables/plan/config/run.yaml` — use that value for Steps 2–4.

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

- All branches: (1) biological context, (2) the full preprocessing plan, (3) **QC review** (after quality filtering and doublet removal, before dimensionality reduction). After QC review, the pipeline runs straight through dimensionality reduction, clustering (fixed resolutions), and UMAP to the final outputs — no further pause.
- On **paired** multiome, the QC review summary also documents the **union doublet removal policy** for confirmation — there is no separate S3 user gate.
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

- **RNA input** (if relevant): full path to `.h5`, MEX directory, or `.h5ad`. Optional **`rna_raw_path`**: companion raw-droplet matrix; when both filtered and raw are supplied, the plan defaults to SoupX for S1a (DecontX when filtered only). Ambient correction is always recommended; the user can opt out at plan review. Optional **`s1a_ambient_method`** in run.yaml (`auto`|`none`|`decontx`|`soupx`).
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

**Execution environment** — if you haven't said already: should this run **locally** on this machine, or on an **HPC cluster** (PBS Pro or SLURM)? For cluster runs I'll probe available queues/partitions on the login node, suggest a project code or account name where I can detect one, and write `deliverables/plan/config/hpc.env` with the `PMA_*` settings you confirm.

Which paths, biological context (if any), and local vs HPC?"

### AGENT_ACTIONS

Once the user answers, in order:

1. Build a minimal `run.yaml` in memory with the fields they supplied:
   ```yaml
   run_dir: <their run dir>
   rna_path: <optional>
   atac_fragments_path: <optional>
   genome_assembly: <required>
   seed: 42
   ```
   Write it to a temporary path like `<run_dir>/run.yaml.draft` (before init copies it to the canonical location).

2. Invoke `executor init --config <draft-run.yaml>`.
   - This scaffolds `internal/` and `deliverables/plan/` (only the plan subtree at init; `figures/`, `qc_review/`, and `results/` appear when outputs are written) and writes the blank biological-context template at `deliverables/plan/config/biological_context.md`.
   - From this point forward, use `$CFG = <run_dir>/deliverables/plan/config/run.yaml` for every CLI call.

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

4. **Configure execution mode + device — confirm with the user first; never auto-default.**
   This is a mandatory one-time gate: `executor run`/`submit` hard-refuse to launch
   any compute until the user's choice is recorded with `--confirmed-by-user`. Ask
   even when local seems obvious — do not assume local just because you're on this
   machine. The gate carries **where to run** (`--mode local|pbs|slurm`) and, for
   **HPC only**, *what device* to configure for the cluster (`--device cpu|gpu`,
   default cpu). **GPU is cluster-only** — never offer `--device gpu` with
   `--mode local`. Preprocessing is CPU-only; `--device gpu` prepares GPU routing for
   the **integration subagent** (future). Explore resources and ask the user — full
   procedure in `inputs_intake.md` Section 4. Only offer `--device gpu` on HPC when a
   cluster GPU partition is detected via `hpc-info`.
   - If the user explicitly chose **local**: `executor configure-execution --config $CFG --mode local --confirmed-by-user`. **Do not run this until the user has actually said local** — there is no "no preference → local" shortcut.
   - If the user said **HPC** (or you're on a login node with `qsub`/`sbatch` and the dataset is large):
     a. Run `executor hpc-info` (parse silently) and measure file sizes with `ls -la` on the input paths the user already provided. Read its GPU fields (`slurm.gpu_partitions` / `suggested_gpu_partition` / `suggested_gpu_gres`) to learn whether GPU is available.
     b. Apply the file-size → scale heuristic from `inputs_intake.md` Section 4 to derive a recommended `PMA_RESOURCES_SCALE`. Select `suggested_partition` / `suggested_account` from `hpc-info` as the candidate values.
     c. Present ONE concrete recommendation (partition + account + scale, **plus a device cpu/gpu choice when a GPU was detected**) with a brief rationale and invite confirmation or override. Do not enumerate the full partition list.
     d. Write settings once confirmed: `executor configure-execution --config $CFG --mode pbs|slurm --pbs-queue ... --pbs-project ... --confirmed-by-user` (or `--slurm-partition` / `--slurm-account`); **add `--device gpu --gpu-gres <suggested_gpu_gres> --gpu-image-uri docker://<registry>/muagene-gpu:<tag>` (SLURM; image_uri may instead come from machine.config) or `--device gpu --pbs-gpu-select-extra 'ngpus=1' --gpu-conda-env muagene-gpu` (PBS — deferred, still conda-env) if the user chose GPU** — see Section 4 for the loud-fail prerequisites. `--confirmed-by-user` records that the user approved the mode; without it, `run`/`submit` will refuse to launch. This records `execution.mode` + `compute.device` + `execution.user_confirmed` in `parameters.yaml` and writes `deliverables/plan/config/hpc.env`.
   - Do not invent partition/account names — use `hpc-info` results. If `hpc-info` returns empty lists or a file is unreachable, see fallback rules in `inputs_intake.md` Section 4.

5. Invoke `executor declare-branch <paired|separate|rna_only|atac_only> --config $CFG`.
   - S0 will confirm this assertion against its own pairing detection; if they mismatch, S0 raises with a clear diff — relay it verbatim.

6. **Run planning (P1 → S0 → P2) — execution location depends on configured mode:**
   The planning phase target is **`plan_review_propose`** (auto-inferred by
   `infer_resume_target()`). `plan_review_propose` depends on `s0_ingest_execute`,
   so Snakemake runs P1 → S0 → P2 → gate-arming in a single head-job submission —
   the gate is armed at the end without a separate step. Targeting `s0_ingest_execute`
   directly would stop one rule early and leave the gate unarmed.
   S0 execution location follows the configured mode (the execution model and
   monitoring mechanics are defined once under *Running on HPC* — don't restate them).
   - **HPC mode (`execution.mode` is `pbs` or `slurm`):** source `hpc.env`, submit the planning head-job, then report-and-repoll:
     ```
     source deliverables/plan/config/hpc.env
     executor submit --config $CFG --executor pbs|slurm
     executor hpc-status --config $CFG             # one-shot: report the daemon's snapshot, then yield
     ```
     Target is auto-inferred (`plan_review_propose`). The Execution-MuAgent daemon
     runs S0 on a compute node and is the sole monitor; follow the **report-and-repoll**
     rule under *Running on HPC* (report one-shot `hpc-status`, then re-poll on a
     non-blocking scheduled wakeup until the gate signal — never block or tail logs).
   - **Local mode (`execution.mode` is `local`):** `executor run --config $CFG --target plan_review_propose` (runs P1 → S0 + plan assembly + gate-arming; ~30s on small inputs).
   - **Do not** retry logic errors (pairing ambiguous, declared-vs-detected mismatch, missing index). Relay the message; user fixes inputs or branch declaration.
   - After S0 completes, surface `validation_report.json` pairing fields if `paired` was downgraded or ambiguous.

### WHAT_TO_SURFACE_BACK

- Confirm `executor init` wrote `deliverables/plan/config/run.yaml` and `biological_context.md`.
- If context was supplied, confirm it was written (`deliverables/plan/config/biological_context.md`).
- If HPC mode was configured, confirm `execution.mode`, the `compute.device` (cpu/gpu) recorded, and the path to `deliverables/plan/config/hpc.env`; remind the user to `source` it before cluster submit/resume. When `device=gpu` on SLURM, also confirm GPU partition/gres and the pinned `gpu_image_uri` (container image pulled from registry — not a conda env).
- After the planning phase completes (`plan_review_propose`), surface `deliverables/plan/context_summary.md` if populated (conflicts or inferred values). Do not paraphrase — paste the markdown back.

See [`stage_prompts/inputs_intake.md`](stage_prompts/inputs_intake.md) for the canonical Step 2 script and the per-context-form handling details.

---

## Step 3 — Confirm the plan

### PROMPT_TO_USER

"Here's the preprocessing plan. Nothing heavy runs until you approve it — this is the single point where cross-stage inconsistencies are cheapest to catch. Review and tell me: **approve as-is**, **revise one or more parameters**, or **abort**.

[... plan_review.md content relayed verbatim ...]"

### AGENT_ACTIONS

1. Run `executor plan-review --intro-context --config $CFG`.
   This prints a JSON object with sample metadata, cell counts, and barcode
   matching data. Do not write anything yet.

2. Write a 100–150-word introductory paragraph using the context data. Rules:
   - Cover all of: organism, tissue, platform/assay type, the aim of the
     analysis (QC → doublet removal → dimensionality reduction → clustering),
     raw cell counts per modality, and the barcode matching result.
   - Write as smooth, user-friendly prose. No bullet points, no jargon.
   - Do not name pipeline stages, step codes, or internal file names
     (e.g., no "S0_ingest", "P2", "preprocessing_plan.json").
   - Use the data exactly as provided; do not round or omit numbers.
   - **Dataset compatibility check (paired multiome candidate only):** This
     check applies only when `workflow_branch = "paired"`. If the barcode
     check found no direct or subset match between RNA and ATAC barcodes
     (i.e., `pairing_confidence` is not "high", or `pairing_status` is not
     "paired"), do not simply flag and ask — instead actively investigate
     and diagnose:
     a. Examine `pairing_ladder` to see which matching steps were attempted
        and why each failed (no overlap, wrong prefix/suffix, ambiguous range).
     b. Cross-check `rna_filtered_status` and `atac_barcodes_source` — a raw
        (non-cell-barcode-filtered) ATAC fragment file will have far more barcodes
        than the RNA cell set, explaining near-zero overlap; check
        `rna_raw_n_barcodes` vs `rna_n_cells` to detect this pattern.
     c. Read `deliverables/plan/config/run.yaml` to inspect the actual
        file paths the user supplied and look for signs of mismatched sources
        (e.g., a filtered matrix paired with unfiltered outputs, samples from a different experiment, or need barcode translations).
     d. Based on (a)–(c), propose the most likely root cause and concrete
        corrective steps (e.g., "re-run with the filtered matrix at
        `filtered_feature_bc_matrix/`, convert bardoes based on direct rules, or supply a `barcode_translation_path` if the two modalities were processed with different whitelists").
     Include this diagnosis and suggestions in the intro paragraph rather
     than a generic warning. Do not apply this check for `rna_only`,
     `atac_only`, or `separate` branches.

3. Invoke `executor plan-review --intro "<paragraph>" --config $CFG`.
   - This re-renders (and writes) BOTH `plan_review_<run>.md` and `plan_summary_<run>.html`
     (the two plan-review deliverables listed in `system_prompt.md`) with the intro
     paragraph prepended before the Summary section.
   - The intro is **persisted**, so the `plan_review_propose` rule and any later
     re-render reproduce it instead of dropping it — you only pass `--intro` once.
   - Also writes per-stage job spec YAMLs to `internal/stage_meta/`.
   - The same content also lives at those paths if the user wants to open them directly.

4. **Marker gene check — mandatory question whenever ambient correction is planned.**
   This is **not** optional agent housekeeping and **not** buried in the approve
   path: surface it as part of presenting the plan. If the plan keeps ambient RNA
   correction (`s1a_ambient.method != none`) **and** the rendered "Marker gene
   expression check" item is still `not set`, you **must** ask the user before any
   approval:
   "The plan runs ambient RNA correction. I recommend checking marker-gene
   expression *before vs after* correction — if a marker shows low ubiquitous
   expression in cells that shouldn't express it, that's ambient contamination, and
   correction should sharpen it back to the right populations. Please give me 5–10
   marker gene symbols to visualise, or tell me to **defer** this to the QC review
   step, or to **skip** it."
   - Escalate the wording to **strongly recommended** when contamination is elevated
     (high `qc_explore` median rho).
   - **Never invent, suggest, look up, or supply gene names yourself** (hard rule,
     [`stage_prompts/qc_threshold_revision.md`](stage_prompts/qc_threshold_revision.md)).
   - The user must make one explicit choice; record it:
     - **provide genes** → `executor revise s1a_ambient "marker_genes=[gene1, gene2, ...]" --config $CFG --rationale "Marker genes provided at plan review"` (stored in `parameters.yaml` as `s1a_ambient.marker_genes`, plotted automatically during S1a).
     - **defer to QC review** → carried as `--defer-marker-genes` on the approve call below (or `--marker-genes defer` on submit).
     - **decline** → carried as `--skip-marker-genes` on the approve call below (or `--marker-genes skip` on submit).
   - If `s1a_ambient.method == none`, skip this question entirely.

5. **QC threshold confirmation — mandatory question, runs after marker-gene step.**
   The "QC strategy" summary item shows `[? needs confirmation]`. Surface it as part of presenting the plan — this is not optional. After the marker-gene question is resolved, ask:
   > "The QC strategy is marked **[? needs confirmation]**. The default MAD-based thresholds are shown in the plan appendix histograms. Would you like to:
   > - Keep the defaults and proceed to approval?
   > - Adjust one or more thresholds (tell me which and how)?
   > - Set an exact threshold value (e.g. RNA `n_genes` lower bound = 300)?
   > - Skip one or more metrics entirely?
   >
   > Any RNA metric (`total_counts`, `n_genes`, `pct_counts_mt`, `pct_counts_ribo`) or ATAC metric (`n_fragments`, `tss_enrichment`, `nucleosome_signal`, `frip`) can be tuned, pinned to an exact value, or disabled — see `qc_threshold_revision.md` → 'Pinning a bound to an exact value' and 'Skipping individual QC metrics'."

   User outcomes:
   - **Confirm defaults** → proceed to approval (step 6 below).
   - **Adjust / skip** → run `executor revise s1_rna_qc <param>=<value> --config $CFG --rationale "<reason>"` (or `s2_atac_qc`). `revise` auto-regenerates the plan deliverables with the revised thresholds and recomputes projected cell-removal counts; re-surface and re-ask until the user confirms.

6. On user decision:
   - **Approve** → `executor approve plan_review --config $CFG --note "approved after review"`, adding `--defer-marker-genes` or `--skip-marker-genes` to match the user's marker-gene choice when no genes were provided. **The executor refuses to approve while the marker-gene decision is unresolved** — if you see that error, you skipped the mandatory question above; go ask it. (On HPC, the same decision is carried as `--marker-genes defer|skip` on `submit --auto-approve`.)
   - **Revise inputs or parameters** → `executor revise <stage> <param>=<value> --config $CFG --rationale "<user's reason>"` (the stage prefix is auto-added to the key, so `revise s1_rna_qc min_counts_floor=500` stores `s1_rna_qc.min_counts_floor=500`). While `plan_review` is unapproved, `revise` **automatically regenerates** `plan_review_<run>.md` / `plan_summary_<run>.html` (and re-derives the cheap QC preview from persisted metrics) so the overlay shows the new value — you do not refresh them by hand. The override wins over the frozen plan when the stage runs. Stage is re-set to awaiting_approval; re-surface the updated deliverable and ask if more revisions are needed before re-approving. (Revise the input knobs — or pin a MAD-derived bound to an exact value with its `*_override` key (e.g. `n_genes_min_override=300`); revising the MAD-derived output key directly has no effect — see [`stage_prompts/qc_threshold_revision.md`](stage_prompts/qc_threshold_revision.md) "Common revise keys" and "Pinning a bound to an exact value".)
   - **Abort** → stop. Tell the user the run dir is intact; they can resume later by re-invoking you on the same config.

### WHAT_TO_SURFACE_BACK

The **Summary** section of `plan_review.md`, verbatim (the appendix is reference detail — surface decision points and any `?` flags from the summary). Do not paraphrase values. Also point the user at `deliverables/plan/plan_summary_<run>.html` (the download-friendly web version — described in the path list in `system_prompt.md`).

If marker genes were stored at this step, confirm the stored gene list in one line (e.g. "Marker genes `Cd3e Cd20 Epcam` stored — the check will run automatically during ambient correction.").

---

## Step 4 — Run with checkpoints

### PROMPT_TO_USER

(At each mandatory pause, surface the stage's proposal content and ask for approval or revision. Pauses are branch-aware:)

- **`p1_context`** (all branches): biological context extraction + conflict resolution. Already handled in Step 2 flow in most cases, but if the user skipped context in Step 2, P1 will stop here.
- **`plan_review`** (all branches): covered in Step 3 — checkpoint **#1**.
- **`post_qc_review`** (all branches): QC review checkpoint **#2** between doublet removal and S4/S5. Generates QC figures in `deliverables/figures/` and `deliverables/qc_review/qc_review_<run_name>.md` (quality-filter and doublet metrics; on **paired**, includes union doublet policy for confirmation). Point the user at `deliverables/qc_review/` for the reports (figures are embedded; raw plots live in `deliverables/figures/`). They may revise thresholds and re-run affected stages before approving. On `separate` / single-modality branches, no cross-modal doublet policy applies. **Hard rule — close the marker-gene loop here:** if `qc_review_<run>.md` contains the notice **"Marker gene expression check not performed"** (this is the second chance when the check was deferred or skipped at plan review), you **must** relay that notice verbatim and obtain an explicit user decision — provide genes → run `executor marker-gene-check --config $CFG <genes...>` (plots before/after and refreshes the QC report), or explicitly decline — **before** approving QC. Do not auto-approve `post_qc_review` past an unaddressed "strongly recommended" notice. Follow [`stage_prompts/qc_threshold_revision.md`](stage_prompts/qc_threshold_revision.md) for the exact procedure. Never supply gene names yourself.
- **`s7_clustering`** (all branches): **not a user checkpoint.** Leiden clustering runs automatically at fixed per-modality resolutions (RNA = 0.7, ATAC = 0.5). Separate / single-modality: these are the **final** `leiden_rna` / `leiden_atac` labels; paired: diagnostic labels for UMAP only. If the user wants different values, `executor revise s7_clustering s7_clustering.rna_resolution=<x>` at plan review.
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
     c. If the stage has a linked summary under `deliverables/qc_review/` (e.g., `qc_review_<run>.md` for post_qc_review), read that too and surface both.
     d. Based on user decision:
        - Approve → `executor approve <stage> --config $CFG`.
        - Revise → if the current stage is `post_qc_review`, follow [`stage_prompts/qc_threshold_revision.md`](stage_prompts/qc_threshold_revision.md) in full; otherwise `executor revise <stage> <key>=<value> --config $CFG`; re-surface the updated proposal; loop.
     e. Re-invoke the appropriate run command (see HPC section below). Continue until `manifest` completes.

**QC threshold revision and post-QC marker gene check:** see [`stage_prompts/qc_threshold_revision.md`](stage_prompts/qc_threshold_revision.md) for the exact HPC and local procedures, artifact deletion rules, plan-vs-`parameters.yaml` behavior, what to surface back, and the procedure for running a marker gene check at QC review when no genes were provided at planning time. At QC review, `marker-gene-check` plots and refreshes reports in one command; use `--plot-only` to skip report refresh.

2. When `manifest` finishes:
   - Read `deliverables/results/run_manifest.json` and extract `workflow_branch`, `outputs`.
   - Point the user at the review notebook `deliverables/results/review_processed_<run>.ipynb` (load + inspect + re-cluster at a custom resolution), the UMAP figures in `deliverables/figures/`, the preprocessing handoff `run_manifest.json`, and the Integration handoff bundle (`post_qc_manifest.json`, `post_qc_<run>.h5mu`). The QC summary lives at `deliverables/qc_review/qc_review_<run>.md`.

### WHAT_TO_SURFACE_BACK

- At each pause: the full proposal yaml content + any linked summary markdown.
- At **`post_qc_review` after a threshold revision:** relay `deliverables/qc_review/qc_review_<run_name>.md` **verbatim** (not the proposal yaml alone). Mention `deliverables/qc_review/qc_summary_<run_name>.html` for the rendered report.
- At completion: the manifest's `outputs` keys + point the user at `post_qc_manifest.json` / `post_qc_<run>.h5mu` (Integration-MuAgent handoff) and a one-line sign-off ("Run complete. Outputs at `deliverables/results/`. I stop here — integration/annotation is out of scope.").

---

## Running on HPC (PBS Pro or SLURM)

The four-step flow above is unchanged on a cluster; only the execution model differs.

**Execution model.** Heavy `_execute` rules dispatch to scheduler jobs; every
`*_propose` rule and the light planning stages (`p1_context`, `plan_review`) plus
`manifest` are `localrules` that run on the orchestrator host (the head-job in
headless mode). **Execution boundary:** `executor run` is local-only and `executor
submit` is cluster-only — Processing-MuAgent never submits or monitors cluster jobs
itself; it prepares the head-job spec + `site.config` and delegates all cluster
execution to Execution-MuAgent via `submit`. There is no `run --executor pbs|slurm`.

**`submit` mechanics.** `submit` writes `internal/stage_meta/head_job.yaml`, then
starts `Execution-MuAgent execute-spec` as a **background supervision daemon** that
submits the head-job, records the job ID to `execution_manifest.jsonl`, and runs the
watch loop (stall detection + kill-on-hang) for the job's lifetime. It returns within
~90 s, once the job ID is confirmed, and **fails loudly if Execution-MuAgent is
absent**. The daemon survives SSH disconnects on most clusters (on
`KillUserProcesses=yes` sites, run `submit` inside `tmux`/`screen`); it writes
`internal/hpc_monitor/monitor_<timestamp>.log` (symlinked `monitor.log`) and removes
`monitor.pid` on exit — that removal is one of the gate signals you wait for.

**S0 ingest** runs inside the planning head-job: in HPC mode, `submit` (no `--target`,
which infers `plan_review_propose`) pulls `s0_ingest_execute` as a dependency on a
compute node — never the login node (QC exploration needs 100+ GB). In local mode,
`run --target plan_review_propose` runs S0 on this machine. Pairing-detection
conflicts are resolved by reading `validation_report.json` after S0 finishes.

### Intake (Step 2) — ask before P1 runs

Execution-mode intake is identical on HPC and local — it lives in **Step 2
AGENT_ACTIONS #4** above and is not repeated here. In short: **always** ask local
vs HPC if the user hasn't said (never default to local silently), and record the
choice with `executor configure-execution ... --confirmed-by-user`; `run`/`submit`
refuse to launch until that confirmation exists. For HPC, probe with `executor
hpc-info`, derive `PMA_RESOURCES_SCALE` from the file-size heuristic, present ONE
partition+account+scale recommendation, then configure. The full heuristic table
and fallback rules (empty partition lists, unreachable files) are in
`inputs_intake.md` Section 4.

### HPC run phases (after plan review)

| Step | Stages | Executes on | You |
|------|--------|-------------|-----|
| Context | P1 | Login node (localrule) | — |
| S0 ingest (+ plan) | S0 | Cluster head-job via `submit` (HPC) / Login node (local) | report `hpc-status`; wait for gate signal |
| Checkpoint **#1** | plan_review | Login node | Review plan |
| QC | S1a → S3 | Cluster head-job via `submit` | report `hpc-status`; wait for gate signal |
| Checkpoint **#2** | post_qc_review | — | Review QC |
| Finish | s_handoff + S4 → S5 → S6 → S7 (clustering) → S8 → manifest | Cluster head-job via `submit` (`s_handoff` is a localrule) | — |

After plan review approval, `source deliverables/plan/config/hpc.env`, then:

- **QC batch:** `executor submit --config $CFG --executor pbs|slurm --auto-approve --auto-approve-except post_qc_review`
- **After QC approval:** `executor submit --config $CFG --executor pbs|slurm` — runs s_handoff + S4→S8→manifest to completion (target `all`, no further gate)

Each gated phase's head-job target is the **gate-arming `*_propose` localrule** (`post_qc_review_propose` for QC), not the phase's last execute stage. Snakemake pulls every execute stage in the phase in as a dependency and runs the propose localrule last, so a single submission runs the whole phase **and** arms the gate. The gate `<stage>` becomes `awaiting_approval` when the propose localrule runs. The final phase (s_handoff + S4→S8→manifest) has no gate, so it targets `all` and runs straight through to the final results in one submission. You never need to run `propose` by hand to surface a gate.

`monitor.pid` removal means the **daemon has stopped** — it is the signal that this phase's compute is over. It *usually* arrives together with the gate sentinel, but not always: if the head job's workflow finishes but its process lingers (an orchestrator that won't exit), the Execution-MuAgent daemon cancels the leftover job and stops, so `monitor.pid` removal may arrive shortly after the gate is already armed. Treat the two as independent signals — which is exactly what the report-and-repoll rule below already does (it drives the next checkpoint on **either** `monitor.pid` gone **or** a gate `awaiting_approval`).

**Job naming convention.** All cluster jobs are named `pma_{stage}_{run_name}` (head jobs) and `pma_{rule}_{run_name}` (child jobs), where `run_name = basename(run_dir)` and `_execute` is stripped from child rule names. Example for run `whelanC57A`: head job = `pma_head_job_whelanC57A`; S1 child = `pma_s1_rna_qc_whelanC57A`. **Always filter squeue by run name** — never use bare `grep pma_head`, which matches jobs from every concurrent sample. The correct form is `squeue -u $USER | grep "pma_head_job_$(basename <run_dir>)"`.

**Monitoring rule — the daemon monitors; Processing reports and re-polls (report-and-repoll).**
After `submit`, the Execution-MuAgent daemon is the sole monitor and writes only structured state. Processing-MuAgent drives this loop:

1. Run one-shot `hpc-status` and report the status to the user in chat.
2. If the job is still running (daemon alive, no gate), **arm a scheduled wakeup for the seconds `hpc-status` printed on its `Next check:` line** — data-driven, typically ~295s (the daemon's 270s heartbeat + a 25s buffer), read from the snapshot's `next_recheck_after_s`. The wakeup's action is to re-run one-shot `hpc-status`.
3. On wake, re-run `hpc-status`. **Report to the user only when the `State:` fingerprint line changed** (scheduler state, health, gate, or findings count); if it is unchanged, stay silent and re-arm. This keeps long jobs quiet while still catching every transition within one heartbeat.
4. **Stop re-arming and drive the next gate** the moment `monitor.pid` is gone or a gate is `awaiting_approval` — `hpc-status` prints `Gate signal present — drive the next checkpoint now`.

A scheduled wakeup is **NON-BLOCKING**: it yields the turn and the runtime re-invokes Processing-MuAgent later (in Claude Code, via `ScheduleWakeup` / `/loop` at the reported cadence). It is **not** a watch loop and **not** `tail -f`. The prohibition on blocking watch loops and `tail -f | grep` remains in force — re-polling via a scheduled wakeup is the *only* sanctioned self-trigger. Never go silent on a state change: the user's only status surface is Processing-MuAgent's chat report via one-shot `hpc-status`.

### How Claude surfaces things during the HPC flow

Behaviour matches local mode until plan review is approved. Beyond the report-and-repoll
rule above, the HPC-specific surfacing is:

- **When the user runs `submit`**, surface the printed PBS/SLURM job ID, the supervision daemon PID, and the log path, then follow the report-and-repoll rule.
- **On approve/revise**, run the appropriate CLI, then `submit` again.

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

These are written to `deliverables/plan/config/hpc.env` by `configure-execution`.

---

## When things go wrong

- **S0 auto-downgrades `paired → separate`** (not an error — this is the new diagnostics-ladder behaviour). User declared `paired` but none of the ladder rungs (direct overlap, suffix-normalized, `barcode_translation_path`, `cell_metadata_path` with `atac_barcode` column) validated cell-level pairing. The run continues on the `separate` branch. Surface `internal/artifacts/s0_ingest/validation_report.json` verbatim — its `pairing.downgrade_reason` field has the specific reason. Ask the user whether to (a) proceed on `separate`, (b) supply a `barcode_translation_path` and rerun S0, or (c) abort and fix inputs upstream (e.g. rerun `cellranger-arc count` for a combined ARC matrix).
- **S0 raises "declared=... conflicts with detected=..."** — the user's `executor declare-branch` doesn't match what S0 detected, AND the declaration is single-modality (`rna_only`/`atac_only`). Single-modality conflicts still raise hard because they signal a data-hygiene problem. Relay the raised message. Ask the user to either correct the declaration or correct the config (drop the unwanted modality).
- **S0 raises "pairing is ambiguous"** — RNA+ATAC Jaccard overlap is between 30% and 80% after normalization/subset checks, and the user did not declare a branch (or declared one that doesn't resolve the ambiguity). Ask the user: are these paired or separate? Based on answer, run `executor declare-branch <paired|separate>` and re-run; for `paired`, supply `barcode_translation_path` if barcode whitelists differ. Don't auto-pick.
- **S3 raises "paired-branch joint barcode intersection is empty"** — S0 committed `paired` but no cell survived both modalities' QC + doublet removal. This usually means QC thresholds were too aggressive. Surface the message; ask the user to revise S1/S2 thresholds via `executor revise s1_rna_qc ...` or `executor revise s2_atac_qc ...`. If the pairing decision used `pairing.translation_table`, also check that the translation table actually covers the QC-surviving cell set.
- **Phase 1 gate raises "biological_context.md is empty"** — the user didn't give context and didn't opt out. Ask for context OR offer the explicit opt-out: `--no-context` on whichever entry point starts the run (`executor run --config $CFG --target plan_review_propose --no-context` in local mode, or `executor submit --config $CFG --executor pbs|slurm --no-context` on HPC).
- **`run`/`submit` raises "Execution mode is not set" or "was not confirmed by the user"** — you tried to launch compute before confirming local vs HPC (this gate fires on fresh runs and resume sessions alike). Stop and confirm the mode with the user: ask local vs HPC, probe `executor hpc-info` for clusters, then record their explicit choice with `executor configure-execution --config $CFG --mode <local|pbs|slurm> --confirmed-by-user`. Never pass `--confirmed-by-user` without having actually asked. Once recorded, re-run the same command and the pipeline proceeds automatically.
- **`run` raises "execution.mode is 'pbs'/'slurm' but `run` is local-only"** — the run is configured for a cluster; `run` only executes locally. Source `hpc.env` and use `executor submit --config $CFG --executor pbs|slurm` instead.
- **S0 OOMs / is Killed / hits walltime in HPC mode** — S0 already runs as a supervised cluster job, so this is a resource-sizing issue, not a location one. Raise `PMA_RESOURCES_SCALE` via `executor configure-execution --config $CFG --mode pbs|slurm --resources-scale N ...`, then `executor submit --config $CFG --executor pbs|slurm` again (omit `--target`). (Re-config of the *same* mode preserves the existing user confirmation — no `--confirmed-by-user` needed for a resource-only change.) **In local mode**, an S0 OOM means the machine is too small — switch to HPC (`configure-execution --mode slurm|pbs`) and submit. There is no automatic local→cluster retry.
- **A stage execute fails at runtime** — relay the failure (HPC: read it from one-shot `executor hpc-status --config $CFG`, which renders the daemon's structured findings; local: snakemake stderr). Do not retry silently; root-cause first. If the user insists on retry, re-`executor submit --executor pbs|slurm --target <stage>_execute` (HPC) or `executor run --config $CFG --target <stage>_execute` (local only).
- **Execution-MuAgent reports `submit_rejected_policy`** — the scheduler rejected the job as a policy error (invalid partition, account, or walltime over the site limit). One-shot `executor hpc-status --config $CFG` renders the scheduler's exact message from the daemon's structured findings. Tell the user which field to correct: partition/account via `executor configure-execution --mode <scheduler> ...` (rewrites `site.config`), or walltime by reducing `PMA_RESOURCES_SCALE`. Then `executor submit` again.
- **Execution-MuAgent reports an environment-preflight error at submit** — provisioning is owned by Execution-MuAgent; `submit` auto-provisions a missing/stale env (policy=auto) but fails loud rather than degrade. Relay the finding verbatim and the fix: `gpu_image_unavailable` → the GPU registry `image_uri` is missing/unreachable (set/fix it; the image is pulled, never built locally); `lock_stale_vs_yaml` → `workflow/envs/processing.yaml` is newer than the lock — run `Processing-MuAgent regenerate-locks` (needs `pip install '.[dev]'`) and commit; `platform_unsupported` → the CPU env is linux-only and this host isn't linux (use a linux host/container); `provision_failed`/`import_failed` → the create/pull or an import failed (relay the stderr tail). On a brand-new machine that was never bootstrapped, run `Execution-MuAgent init-machine --processing-repo <repo>` first.
- **Per-stage specs not written** — specs are written automatically by `executor plan-review`. If `internal/stage_meta/` is missing or empty, re-run `executor plan-review --config $CFG`. Specs are internal state; do not surface them to the user unless asked.
- **`hpc-status` shows "Supervisor: not running" alongside a RUNNING or PENDING scheduler state** — the supervision daemon has died but the cluster job is still active. Without the daemon, stalled jobs will not be auto-cancelled. Restart the daemon: `executor supervisor-restart --config $CFG`. This resumes the full watch loop (stall detection, kill-on-hang) against the already-running job without resubmitting. Tell the user what happened and what you did.
- **Supervision daemon crashes on a site with KillUserProcesses=yes** — when the user's SSH session ends, systemd kills all their processes including the daemon. The cluster job keeps running, but protection is gone. For the current run, tell them to use `supervisor-restart` as soon as they reconnect. Going forward, suggest running `submit` inside a `tmux` or `screen` session on that cluster.
- **One-shot `hpc-status` shows "review gate awaiting approval" before the pipeline has made any progress** — a stale `awaiting_approval` sentinel from a prior run (or an old head job still writing to it) is blocking progress. It now shows up directly in one-shot status before any real progress. Fix: (1) `squeue -u $USER | grep "pma_head_job_$(basename <run_dir>)"` → `scancel <JOBID>` for each result; (2) `rm internal/proposals/<stage>.awaiting_approval`; (3) re-run `executor submit` and report the next one-shot `hpc-status`.
- **Tempted to monitor a long-running job yourself** — don't. Rely on the daemon (the sole monitor) and read one-shot `executor hpc-status --config $CFG`; never run a blocking loop or `tail -f | grep`. (Re-polling via a non-blocking scheduled wakeup per the report-and-repoll rule is not a blocking loop and is the sanctioned way to re-check.)
- **Blank ATAC QC figures ("(no data)") / `qc_explore` log shows `chrom_bound_filter_failed: bgzip not found on PATH`** — an *execution-environment* error, not a scientific one. The cluster child job did not have the project conda env's tools (`bgzip`/`tabix`, which live in `$PMA_CONDA_ENV/bin`) on PATH, so the ATAC fragment chr-renaming + chromosome-bound filter was skipped, the SnapATAC2 import matched zero fragments (Ensembl-named fragments vs UCSC-named reference), and `atac_qc_metrics.parquet` came back empty. **This is fixed structurally:** each generated Snakemake child jobscript is sanitized to self-activate `$PMA_CONDA_ENV` (`executor.hpc.sanitize_snakemake_jobscript` → `inject_conda_activation_text`, invoked by `slurm-submit.sh` / `pbs-submit.sh`), and `htslib` is pinned in `workflow/envs/processing.yaml` for `--use-conda` runs. **Note the default cluster profiles set `use-conda: false`**, so per-rule conda envs are *not* built — protection then depends entirely on the `--conda-env` you configured actually containing `bgzip`/`tabix`. Reading gzipped fragments uses Python's `gzip` (no `bgzip` needed for reads); the BED4→fragments conversion (`executor/io.py`) *requires* `bgzip`/`tabix` and **raises a clear, actionable error** if they are absent — and an empty ATAC import now raises rather than silently emitting a blank figure (it does **not** silently fall back). If you hit it: confirm `configure-execution --conda-env <name>` was set (recorded in `hpc.env` / `site.config`) and that `bgzip`/`tabix` exist in that env (`conda run -n <env> bgzip --version`).
