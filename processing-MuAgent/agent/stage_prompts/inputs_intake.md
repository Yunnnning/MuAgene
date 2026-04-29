# Step 2 — Inputs intake (paths + optional biological context)

Script for the turn(s) after the user has declared their analysis type in Step 1. Goal: collect enough to build a valid `run.yaml`, write it via `executor init`, populate biological context (if offered), declare the branch, and run up to the plan-review gate.

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
> If the RNA and ATAC barcodes overlap ≥99% I'll treat as `paired`; if they don't overlap at all I'll treat as `separate`. Anything in between, I'll stop and ask.

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
- `<run_dir>/deliverables/{pre_run,post_run}/summary/`, plus `post_run/{figures,processed,notebooks}/`

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

### 4. Declare the branch

```
executor declare-branch <rna_only|atac_only|paired|separate> --config $CFG
```

This writes `plan.workflow_branch_declared` to `parameters.yaml` as a `source=user` assertion. S0 will confirm it against its own pairing detection and raise with a clear diff if they don't match.

### 5. Run to the plan-review gate

```
executor run --config $CFG --target p2_plan_execute
```

This triggers:

- `p1_context_propose` → `p1_context_execute` (extracts context fields, fetches DOIs)
- `s0_ingest_propose` → `s0_ingest_execute` (validates inputs, detects pairing, confirms branch)
- `p2_plan_propose` → `p2_plan_execute` (assembles the full plan)

Stops at the `plan_review` gate. Takes ~30s on the example data; up to a few minutes on larger inputs.

## What to surface back

After `executor init`: confirm the canonical config path and blank context template path.

After biological-context write (cases a/b/c): confirm the file exists at `deliverables/pre_run/config/biological_context.md` and that you populated the fields the user told you about.

After `executor run --target p2_plan_execute`:

- If `deliverables/pre_run/summary/context_summary.md` exists, paste its content back verbatim. Any conflicts (e.g. "report says mouse, file fingerprint says human") surface here and must be resolved before Step 3.
- If `executor run` errored:
  - **Phase 1 gate error** — context template is blank and user didn't opt out. Ask for context OR tell them to re-invoke with `--no-context`.
  - **S0 declared-vs-detected mismatch** — relay the raised error and ask the user to fix the declaration or the inputs.
  - **S0 ambiguous pairing** — relay the raised error; ask paired vs separate; re-run `executor declare-branch` and re-try.

Transition to Step 3 (`plan_review`) once `plan_summary.md` and `plan_review.md` exist under `deliverables/pre_run/summary/`.

## Explicit non-actions

- Do NOT write `parameters.yaml` directly — use `executor declare-branch` / `executor revise`.
- Do NOT touch `biological_context.md` directly — always route through `context_mapper.write_report`.
- Do NOT copy input files into the run dir. Paths are referenced in place.
- Do NOT auto-retry on S0 errors. Relay the message; let the user correct the root cause.
- Do NOT skip biological context "to save time" — the Phase 1 gate exists because context shapes downstream QC thresholds; the cost of silently proceeding is recovering later from wrong defaults.
