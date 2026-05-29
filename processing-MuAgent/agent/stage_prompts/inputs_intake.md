# Step 2 — Inputs intake (paths + optional biological context)

Script for the turn(s) after the user has declared their analysis type in Step 1. Goal: collect enough to build a valid `run.yaml`, write it via `executor init`, populate biological context (if offered), configure execution mode (local vs HPC), declare the branch, and run up to the plan-review gate.

## What to say

Tailor the required paths to the declared `workflow_branch`:

### For `rna_only`

> I need:
> - **RNA input path** — one of: 10x Cell Ranger `.h5`, 10x MEX directory (matrix.mtx + barcodes.tsv + features.tsv), or `.h5ad`.
> - **Genome assembly** — e.g. `mm10`, `GRCh38`. Required; no default. I'll cross-check it against features in the matrix where I can.
>
> Optional:
> - **Biological context** — organism, tissue, assay, any DOIs. Free text is fine, or paste a filled Biological Context Report, or give me a path to one.
> - **Study goal** — `clustering_inference` (default) or `rare_populations`. Shapes my S3 doublet-policy recommendation.
> - **Seed** — default 42.
> - **Execution** — Should I run locally on this machine, or submit jobs to a cluster (HPC: PBS Pro or SLURM)? - If you choose HPC and have not yet set the required `PMA_*` environment variables, I'll run `hpc-info` on the login node, list available queues/partitions, suggest a project code or account where I can detect one, and ask you to confirm before I write `hpc.env`.

### For `atac_only`

> I need:
> - **ATAC fragments path** — `fragments.tsv.gz`. The `.tbi` index must sit next to it with the same stem.
> - **Genome assembly** — e.g. `mm10`, `GRCh38`. Required. I cross-check the fragments' chromosome names against it (`chr1..chrY,chrM` vs `1..22,X,Y,MT`).
>
> Optional: same as rna_only.

### For `paired` or `separate`

> I need:
> - **RNA input path** (h5 / MEX / h5ad)
> - **ATAC fragments path** (`fragments.tsv.gz` + matching `.tbi`)
> - **Genome assembly**
>
> Optional: same as above.
>
> If the RNA and ATAC barcodes have Jaccard overlap ≥80%, or one modality's barcodes are ≥80% contained in the other (typical when cell counts differ before QC), I'll treat as `paired`. If they don't overlap at all I'll treat as `separate`. Jaccard between 30% and 80% with no subset relation, I'll stop and ask.

## Actions

Once the user answers, execute in order:

### 1. Draft `run.yaml`

Build an in-memory dict from the user's answers (omit fields not supplied):

```yaml
run_dir: <user's run_dir>
rna_path: <optional>
atac_fragments_path: <optional>
genome_assembly: <required>
study_goal: <clustering_inference | rare_populations>     # default clustering_inference
seed: 42
```

Write this to a draft path, e.g. `<run_dir>/run.yaml.draft`, or to any writable location you choose. The draft is read once by `executor init` and copied to the canonical location; its location after init doesn't matter.

### 2. `executor init`

```
executor init --config <draft-run.yaml>
```

This creates:

- `<run_dir>/internal/` — pipeline state scaffold
- `<run_dir>/deliverables/pre_run/config/run.yaml` — canonical copy of the config
- `<run_dir>/deliverables/pre_run/config/biological_context.md` — blank template
- `<run_dir>/deliverables/{pre_run,checkpoint,post_run}/` — user-facing outputs split by lifecycle phase

From now on, `$CFG = <run_dir>/deliverables/pre_run/config/run.yaml` for every subsequent CLI call.

### 3. Populate biological context

Three cases based on what the user gave you:

**(a) Short chat text** — you have free-form text about organism/tissue/assay + maybe a DOI. Do NLU yourself (you're the LLM) and extract structured fields; pass them to `context_mapper`:

```python
from executor import context_mapper
md = context_mapper.build_report_from_chat(
    organism="mouse",                                    # extracted from user text
    tissue="testis",
    assay="single-nucleus multiome (snRNA + snATAC)",
    dois=["10.1016/j.stemcr.2025.102449"],                # any DOIs they mentioned
    notes="GSE268104; adult mouse C57BL/6",               # free-form extras
)
context_mapper.write_report(run_dir, md)
```

Empty strings for fields the user didn't mention — don't invent values. `is_unfilled_template` will only flag the report as unfilled if all three of organism/tissue/assay are empty.

**(b) Filled template path** — user gave `/path/to/their_context.md`. Read the file and pass its content unchanged:

```python
from executor import context_mapper
content = Path("/path/to/their_context.md").read_text()
context_mapper.write_report(run_dir, content)
```

**(c) DOI list only** — user gave just DOIs. Merge them into the blank template:

```python
from executor import context_mapper
md = context_mapper.build_report_from_chat(dois=["10.xxxx/...", "10.yyyy/..."])
context_mapper.write_report(run_dir, md)
```

P1 will fetch abstracts for each DOI during `executor run --target p2_plan_execute`.

**(d) Nothing supplied** — leave the blank template. Warn the user that the Phase 1 gate will block them and they'll need to either paste context later or opt out explicitly with `executor run --config $CFG --no-context`. Don't opt out silently on their behalf.

### 4. Configure execution mode

If the user has not already said **local** vs **HPC**, ask now (before P1 runs).

**Local (default for small data / laptop):**
```
executor configure-execution --config $CFG --mode local
```

**HPC (PBS or SLURM):**

1. Run `executor hpc-info` and paste the JSON back to the user in a readable summary:
   - `detected_scheduler` (`pbs` | `slurm`)
   - available `pbs.queues` or `slurm.partitions`
   - `pbs.suggested_project` / `slurm.suggested_account` (if any)
   - `current_env` — any `PMA_*` vars already exported
2. Ask the user to confirm or override:
   - PBS: queue (+ project code if their site requires `-P`)
   - SLURM: partition (+ account if required)
   - `PMA_NOTIFY_EMAIL` (recommended)
   - optional `PMA_RESOURCES_SCALE`
3. Write settings once confirmed:
   ```
   executor configure-execution --config $CFG --mode pbs \
       --pbs-queue <queue> --pbs-project <project> --notify-email <email>
   ```
   (or `--mode slurm --slurm-partition ... --slurm-account ...`)

Do **not** guess queue/partition/project/account — use `hpc-info` suggestions and wait for user confirmation. If `hpc-info` returns empty lists, ask the user for the values directly.

### 5. Declare the branch

```
executor declare-branch <rna_only|atac_only|paired|separate> --config $CFG
```

This writes `plan.workflow_branch_declared` to `parameters.yaml` as a `source=user` assertion. S0 will confirm it against its own pairing detection and raise with a clear diff if they don't match.

### 6. Run to the plan-review gate (flexible S0)

**Default — try local first:**

```
executor run --config $CFG --target p2_plan_execute
```

Runs P1 → S0 → P2 and stops at `plan_review`. Small inputs: ~30s.

**Large dataset upfront** (user says so, or dense_txt / very large h5): configure HPC first, then:

```
source deliverables/pre_run/config/hpc.env
executor run --config $CFG --executor pbs|slurm --target s0_ingest_execute
executor run --config $CFG --target p2_plan_execute
```

**S0 failed locally with a resource error** (OOM, Killed, walltime — check with
`from executor.hpc import looks_like_resource_failure` on snakemake stderr):

1. Configure HPC if not done (`hpc-info` + `configure-execution`; raise `PMA_RESOURCES_SCALE` if needed).
2. `source deliverables/pre_run/config/hpc.env`
3. `executor run --config $CFG --executor pbs|slurm --target s0_ingest_execute`
4. `executor run --config $CFG --target p2_plan_execute`

Do **not** cluster-retry logic errors (pairing ambiguous, path missing, branch mismatch). Relay and let the user fix inputs or `declare-branch`.

## What to surface back

After `executor init`: confirm the canonical config path and blank context template path.

After biological-context write (cases a/b/c): confirm the file exists at `deliverables/pre_run/config/biological_context.md` and that you populated the fields the user told you about.

After `configure-execution`: confirm `execution.mode` and, for HPC, the path to `deliverables/pre_run/config/hpc.env`. Tell the user to `source` that file before any cluster submit/resume.

After `executor run --target p2_plan_execute`:

- If `deliverables/pre_run/summary/context_summary.md` exists, paste its content back verbatim. Any conflicts (e.g. "report says mouse, file fingerprint says human") surface here and must be resolved before Step 3.
- If `executor run` errored:
  - **Phase 1 gate error** — context template is blank and user didn't opt out. Ask for context OR tell them to re-invoke with `--no-context`.
  - **S0 declared-vs-detected mismatch** — relay the raised error and ask the user to fix the declaration or the inputs.
  - **S0 ambiguous pairing** — relay the raised error; ask paired vs separate; re-run `executor declare-branch` and re-try.

Transition to Step 3 (`plan_review`) once `plan_review.md` exists under `deliverables/pre_run/summary/` (written by `plan_review_propose` or `processing-muagent plan-review`).

## Explicit non-actions

- Do NOT write `parameters.yaml` directly — use `executor declare-branch` / `executor revise`.
- Do NOT touch `biological_context.md` directly — always route through `context_mapper.write_report`.
- Do NOT copy input files into the run dir. Paths are referenced in place.
- Do NOT auto-retry on S0 errors. Relay the message; let the user correct the root cause.
- Do NOT skip biological context "to save time" — the Phase 1 gate exists because context shapes downstream QC thresholds; the cost of silently proceeding is recovering later from wrong defaults.
