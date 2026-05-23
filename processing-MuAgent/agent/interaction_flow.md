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

For any of these, I do QC → doublet detection → dim reduction → clustering → UMAP per modality, and then I stop. I don't do integration, annotation, or any downstream analysis.

Before you answer, one other thing I'll need: a run directory — a writable folder where I can put intermediate artifacts, approval checkpoints, and final outputs.

Mandatory approval points depend on the branch:

- All branches: (1) biological context, (2) the full preprocessing plan, (3) clustering resolution.
- `paired` / `separate` only: (4) doublet-policy reconciliation at S3 — how to combine the RNA and ATAC detector calls.
- `rna_only` / `atac_only`: S3 runs with the single-detector recommended policy; no reconciliation to confirm, so no pause there unless you ask for stage-by-stage review.

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

- **RNA input** (if relevant): full path to `.h5`, MEX directory, or `.h5ad`.
- **ATAC input** (if relevant): full path to `fragments.tsv.gz`. The `.tbi` must sit next to it — I'll fail fast if it doesn't.
- **Genome assembly**: `mm10`, `GRCh38`, etc. Required. I cross-check this against the ATAC fragment chromosome naming, so declare it carefully.
- **Seed** (optional, default 42).

Biological context is optional but strongly recommended — it shapes QC thresholds and surfaces conflicts early. You can give it to me in any of three forms:

1. **Short chat text** — just tell me: organism, tissue, assay, any related DOIs.
2. **A filled Biological Context Report** — paste the content, or give me a path to an existing markdown file.
3. **A DOI list only** — I'll fetch abstracts and extract what I can.

Which paths, and do you want to include biological context?"

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
   - This scaffolds `internal/` and `deliverables/{pre_run,post_run}/` and writes the blank biological-context template at `deliverables/pre_run/config/biological_context.md`.
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

4. Invoke `executor declare-branch <paired|separate|rna_only|atac_only> --config $CFG`.
   - S0 will confirm this assertion against its own pairing detection; if they mismatch, S0 raises with a clear diff — relay it verbatim.

5. Invoke `executor run --config $CFG --target p2_plan_execute`.
   - This runs P1 context extraction + S0 validation + P2 plan assembly, stopping at the `plan_review` gate. Takes ~30s on small inputs.

### WHAT_TO_SURFACE_BACK

- Confirm `executor init` wrote `deliverables/pre_run/config/run.yaml` and `biological_context.md`.
- If context was supplied, confirm it was written (`deliverables/pre_run/config/biological_context.md`).
- After `executor run --target p2_plan_execute`, surface `deliverables/pre_run/summary/context_summary.md` if populated (conflicts or inferred values). Do not paraphrase — paste the markdown back.

See [`stage_prompts/inputs_intake.md`](stage_prompts/inputs_intake.md) for the canonical Step 2 script and the per-context-form handling details.

---

## Step 3 — Confirm the plan

### PROMPT_TO_USER

"Here's the preprocessing plan. Nothing heavy runs until you approve it — this is the single point where cross-stage inconsistencies are cheapest to catch. Review and tell me: **approve as-is**, **revise one or more parameters**, or **abort**.

[... plan_review.md content relayed verbatim ...]"

### AGENT_ACTIONS

1. Invoke `executor plan-review --config $CFG`.
   - This re-renders (and writes) `deliverables/pre_run/summary/plan_review.md` — the concise 8-item review summary.
   - The same content also lives at that path if the user wants to open it directly.

2. On user decision:
   - **Approve** → `executor approve plan_review --config $CFG --note "approved after review"`.
   - **Revise a parameter** → `executor revise <stage> <key>=<value> --config $CFG --rationale "<user's reason>"`. Stage is re-set to awaiting_approval; ask if more revisions are needed before re-approving.
   - **Abort** → stop. Tell the user the run dir is intact; they can resume later by re-invoking you on the same config.

### WHAT_TO_SURFACE_BACK

The full 8-item `plan_review.md` content, verbatim. Do not paraphrase values. If the review flags `?` uncertainties or unresolved warnings, highlight those — they're the decision points.

---

## Step 4 — Run with checkpoints

### PROMPT_TO_USER

(At each mandatory pause, surface the stage's proposal content and ask for approval or revision. Pauses are branch-aware:)

- **`p1_context`** (all branches): biological context extraction + conflict resolution. Already handled in Step 2 flow in most cases, but if the user skipped context in Step 2, P1 will stop here.
- **`plan_review`** (all branches): covered in Step 3.
- **`s3_doublets`** (`paired` / `separate` only): Scrublet + ATAC detector overlap table; user confirms reconciliation policy (union vs intersection). For `rna_only` / `atac_only`, auto-approve silently — no reconciliation to confirm — unless the user explicitly asked for per-stage review.
- **`s7_clustering`** (all branches): resolution sweep results; user confirms per-modality resolution or revises.

(For other stages, auto-approve silently unless the user asked for per-stage review.)

"I'm at the **[stage]** checkpoint. Here's what the pipeline proposes:

[... proposal yaml + summary markdown content, verbatim ...]

Approve, revise, or abort?"

### AGENT_ACTIONS

1. Invoke `executor run --config $CFG` (no `--auto-approve`).
   - Snakemake runs every stage whose approval sentinel exists, stopping at the first stage whose `.approved` is missing.
   - Loop:
     a. Run `executor status --config $CFG` to see which stage is currently `awaiting_approval`.
     b. Read `<run_dir>/internal/proposals/<stage>.yaml` — the structured proposal.
     c. If the stage has a linked summary in `deliverables/post_run/summary/` (e.g., `resolution_summary.md` for s7), read that too and surface both.
     d. Based on user decision:
        - Approve → `executor approve <stage> --config $CFG`.
        - Revise → `executor revise <stage> <key>=<value> --config $CFG`; re-surface the updated proposal; loop.
     e. Re-invoke `executor run --config $CFG`. Continue until `manifest` completes.

2. When `manifest` finishes:
   - Read `deliverables/post_run/summary/run_manifest.json` and extract `workflow_branch`, `outputs`.
   - Point the user at `deliverables/post_run/summary/qc_summary.md`, the final figures in `deliverables/post_run/figures/`, and the handoff artifact `run_manifest.json`.

### WHAT_TO_SURFACE_BACK

- At each pause: the full proposal yaml content + any linked summary markdown.
- At completion: the manifest's `outputs` keys + a one-line sign-off ("Run complete. Outputs at `deliverables/post_run/`. I stop here — integration/annotation is out of scope.").

---

## Running on HPC (PBS Pro or SLURM)

The four-step flow above is unchanged on a cluster. Only the underlying execution
model differs: heavy `_execute` rules dispatch to scheduler jobs, while every
`*_propose` rule, the planning stages (`p1_context`, `p2_plan`, `plan_review`),
`s0_ingest`, and `manifest` are declared `localrules` and run on the orchestrator
host (login node in interactive mode; the head-job in headless mode).

### Three phases

1. **Planning (Phase A)** — on the login node, inside a `tmux`/`screen` session:
   Steps 1–3 above plus `s0_ingest`. S0 is local because it can raise a
   pairing-detection conflict that needs interactive resolution before any
   cluster job is dispatched. The user reviews `pre_run` deliverables here
   (`plan_review.md`, `context_summary.md`).

   ```bash
   # Inside tmux on the login node:
   processing-muagent init --config config/run.yaml
   processing-muagent run --config $CFG --target s0_ingest_execute
   ```

2. **Main preprocessing (Phase B)** — submit the unattended head-job.
   `s3_doublets` defaults to the union policy locked in at `plan_review`; the
   head-job runs S1a → S6 → S7_propose then stops because
   `s7_clustering.approved` is missing. The exit hook emails `$PMA_NOTIFY_EMAIL`
   if set.

   ```bash
   processing-muagent submit --config $CFG --executor pbs \
       --auto-approve --auto-approve-except s7_clustering
   ```

3. **Resolution review + finish (Phase C)** — open the static review HTML at
   `<run_dir>/deliverables/post_run/notebooks/resolution_review.html`
   (no Jupyter needed). Approve or revise, then submit a small finishing job
   that runs `s7_clustering_execute` → `s8_umap` → `manifest`.

   ```bash
   # Approve:
   processing-muagent approve s7_clustering --config $CFG
   # OR revise:
   processing-muagent revise s7_clustering s7_clustering.rna.resolution=1.2 --config $CFG

   processing-muagent submit --config $CFG --executor pbs
   ```

### How Claude surfaces things during the HPC flow

- **At Step 1–4 of the interactive flow**, behaviour is identical to local mode.
- **When the user runs `submit`**, surface the printed PBS/SLURM job id and remind
  them they can poll with `processing-muagent status --watch --config $CFG`.
- **When the email arrives (or `status --watch` shows `s7_clustering
  awaiting_approval`)**, surface the path to `resolution_review.html` (the
  primary review artifact) AND `resolution_review.ipynb` (for power users who
  want to re-cluster at custom resolutions interactively). Paste the contents
  of `resolution_summary.md` verbatim, plus the proposal's `review_artifacts`
  block.
- **On approve/revise**, run the appropriate CLI as today, then `submit` again
  (or `run --executor pbs` if the user prefers foreground in a tmux).

### Site setup (one time)

```bash
# Imperial RDS (PBS Pro):
export PMA_PBS_QUEUE=v1_throughput72         # or your allocation's queue
export PMA_PBS_PROJECT=<your project code>   # if your queue needs -P
export PMA_NOTIFY_EMAIL=<you@example.com>

# Generic SLURM site:
export PMA_SLURM_PARTITION=cpu
export PMA_SLURM_ACCOUNT=<your account>
export PMA_NOTIFY_EMAIL=<you@example.com>

# Optional: scale memory + walltime for large datasets (default 1).
export PMA_RESOURCES_SCALE=2                  # ~30k-cell, =4 for ~100k-cell
```

---

## When things go wrong

- **S0 raises "declared=... conflicts with detected=..."** — the user's `executor declare-branch` doesn't match what S0 detected. Relay the raised message. Ask the user to either correct the declaration (re-run Step 2 from `executor declare-branch`) or correct the config (edit paths and re-run from Step 2).
- **S0 raises "pairing is ambiguous"** — RNA+ATAC barcode overlap is between 30% and 99% after normalization. Ask the user: are these paired or separate? Based on answer, run `executor declare-branch <paired|separate>` and re-run. Don't auto-pick.
- **Phase 1 gate raises "biological_context.md is empty"** — the user didn't give context and didn't opt out. Ask for context OR offer `executor run --config $CFG --no-context` as the explicit opt-out.
- **A stage execute fails at runtime** — relay the traceback from snakemake. Do not retry silently; root-cause first. If the user insists on retry, use `executor run --config $CFG --target <stage>_execute`.
