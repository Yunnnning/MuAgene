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
study_goal: <clustering_inference | rare_populations>     # default clustering_inference; guides S1a + S3
# s1a_ambient_method: <auto | none | decontx | soupx>   # optional; else plan uses study_goal + ingest
seed: 42
```

Write this to a draft path, e.g. `<run_dir>/run.yaml.draft`, or to any writable location you choose. The draft is read once by `executor init` and copied to the canonical location; its location after init doesn't matter.

### 2. `executor init`

```
executor init --config <draft-run.yaml>
```

This creates:

- `<run_dir>/internal/` — pipeline state scaffold
- `<run_dir>/deliverables/plan/config/run.yaml` — canonical copy of the config
- `<run_dir>/deliverables/plan/config/biological_context.md` — blank template
- `<run_dir>/deliverables/plan/` — created at init; `figures/`, `qc_review/`, and `results/` appear when outputs are written

From now on, `$CFG = <run_dir>/deliverables/plan/config/run.yaml` for every subsequent CLI call.

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

P1 will fetch abstracts for each DOI during the planning phase (`executor run --target plan_review_propose` locally, or `executor submit` on HPC — both pull P1 → S0 and assemble the plan in-process).

**(d) Nothing supplied** — leave the blank template. Warn the user that the Phase 1 gate will block them and they'll need to either paste context later or opt out explicitly with `executor run --config $CFG --no-context`. Don't opt out silently on their behalf.

### 4. Configure execution mode

If the user has not already said **local** vs **HPC**, ask now (before P1 runs) —
**always confirm; never auto-default to local.** `executor run`/`submit` hard-refuse
to launch any compute until the user's choice is recorded with `--confirmed-by-user`
(this gate fires on fresh runs and resume sessions alike). It is a one-time gate:
once confirmed, the rest of the pipeline runs automatically.

**This gate carries two choices — explore the resources, then ask the user:**
- **Where to run** — `--mode local | pbs | slurm` (always ask; never auto-default).
- **What device (HPC only, integration subagent)** — `--device cpu | gpu` (default `cpu`).
  Preprocessing stages are **CPU-only** (`_GPU_CAPABLE` is empty). `--device gpu` on HPC
  prepares cluster GPU infrastructure (container pull, partition/gres routing) for the
  **integration subagent** to use later — it does not accelerate preprocessing today.

**Local** (only after the user explicitly chooses it — do not assume local just
because you're on this machine):
```
executor configure-execution --config $CFG --mode local --confirmed-by-user
```
Do not pass `--device gpu` with local mode — `configure-execution` rejects it.

**HPC (PBS or SLURM):**

1. Run `executor hpc-info`. Parse the JSON silently — do not dump the raw JSON to the user.
   Also read the **GPU** fields: `slurm.gpu_partitions`, `slurm.suggested_gpu_partition`,
   `slurm.suggested_gpu_gres` (SLURM). A non-empty `gpu_partitions` / a `suggested_gpu_gres`
   means GPU is available on this cluster — surface it as a device choice in step 4. (PBS GPU
   syntax is site-variable, so `hpc-info` does not auto-suggest a PBS GPU select; ask the user
   for the `ngpus=…` form if they want GPU on PBS.)

2. Measure input file sizes. For every path the user already provided this turn, run:
   ```bash
   ls -la <path>
   ```
   If the path is a directory (MEX format), measure `matrix.mtx.gz` inside it:
   ```bash
   ls -la <rna_dir>/matrix.mtx.gz
   ```
   If a path is unreachable (permission error, NFS timeout, file not found), note it and apply the fallback rules below.

3. Apply the scale heuristic silently:

   **RNA inputs** (`.h5`, `matrix.mtx.gz`, `.h5ad`):
   | File size        | Estimated cells | Recommended `PMA_RESOURCES_SCALE` |
   |------------------|-----------------|-----------------------------------|
   | < 100 MB         | ~1–10 k         | 1                                 |
   | 100 MB – 500 MB  | ~10–50 k        | 2                                 |
   | > 500 MB         | ~50 k+          | 4                                 |

   **ATAC inputs** (`fragments.tsv.gz`):
   | File size       | Relative size | Recommended `PMA_RESOURCES_SCALE` |
   |-----------------|---------------|-----------------------------------|
   | < 300 MB        | small         | 1                                 |
   | 300 MB – 1 GB   | medium        | 2                                 |
   | > 1 GB          | large         | 4                                 |

   When both modalities are present, take the **maximum** of the two recommended scales.

   **Fallback rules:**
   - File unreachable → note it in the recommendation, ask the user to confirm or supply scale manually.
   - `hpc-info` returns empty `pbs.queues` / `slurm.partitions` → note "no partitions detected", ask the user to supply the value directly.
   - `hpc-info` returns no `suggested_account` / `suggested_project` → omit from recommendation; ask if the site requires one.

4. Present ONE concrete recommendation to the user. Include a **Device** line whenever
   `hpc-info` showed GPU availability — present cpu vs gpu as an explicit choice; do not
   default to gpu. Example format:
   > Based on your RNA input (~180 MB, ~10–50 k cells), I recommend:
   > - **Partition:** cpu (detected from your cluster)
   > - **Account:** project_abc (detected from environment)
   > - **Scale:** `PMA_RESOURCES_SCALE=2`
   > - **Device:** your cluster has a GPU partition (`gpu`, `gpu:A5000:1`). GPU is for
   >   the **integration subagent** (future) — preprocessing stays on CPU. Do you want
   >   to configure GPU routing now (`--device gpu`), or keep everything on **CPU**?
   >
   > Does this look right, or would you like to change any of these?

   Adapt the wording to what was actually detected. Omit the Device line entirely when no GPU
   was detected. Do not enumerate the full partition or account list unless the user asks —
   just state the chosen values with a one-line rationale.

5. Write settings once the user confirms (or overrides):
   ```
   executor configure-execution --config $CFG --mode pbs \
       --pbs-queue <queue> --pbs-project <project> --confirmed-by-user
   ```
   (or `--mode slurm --slurm-partition ... --slurm-account ... --confirmed-by-user`).
   `--confirmed-by-user` records the user's approval; without it `run`/`submit` refuse to launch.

   **If the user chose GPU**, add the device flags (sourced from `hpc-info`'s GPU fields):
   ```
   executor configure-execution --config $CFG --mode slurm \
       --slurm-partition <cpu_partition> --slurm-account <account> \
       --device gpu --gpu-partition <suggested_gpu_partition> --gpu-gres <suggested_gpu_gres> \
       --gpu-image-uri docker://<registry>/muagene-gpu:<tag> --confirmed-by-user
   ```
   (PBS GPU is deferred — still uses a GPU conda env: `--device gpu --pbs-gpu-select-extra 'ngpus=1' [--gpu-queue <q>] --gpu-conda-env muagene-gpu`.)
   `configure-execution` fails loud on missing prerequisites — pre-empt them: SLURM `--device gpu`
   **requires** `--gpu-gres` and `--gpu-image-uri` (the SLURM GPU env is a container PULLED from that
   pinned reference — or set `gpu_image_uri` once in `~/.muagene/machine.config` via init-machine);
   PBS defaults to `ngpus=1` if you omit the select-extra.
   Add `--singularity-module <module>` when the site needs `module load` for singularity.

Do **not** invent partition/account names — use `hpc-info` results only. If `hpc-info` returns empty lists, ask the user for the values directly.

### 5. Declare the branch

```
executor declare-branch <rna_only|atac_only|paired|separate> --config $CFG
```

This writes `plan.workflow_branch_declared` to `parameters.yaml` as a `source=user` assertion. S0 will confirm it against its own pairing detection and raise with a clear diff if they don't match.

### 6. Run to the plan-review gate

The planning target is **`plan_review_propose`** (auto-inferred by `submit` when
`--target` is omitted). `plan_review_propose` depends on `s0_ingest_execute`, so
Snakemake runs P1 → S0 → plan assembly → gate-arming in one invocation — the
`plan_review` gate is armed at the end without a separate step. Do **not** target
`s0_ingest_execute` alone: that stops one rule early and leaves the gate unarmed.

S0 (merged planning compute: load + validate + pair + QC explore) runs inside that
dependency chain. The former separate `p2_plan` rule was merged into S0 — there is no
`p2_plan_execute` rule; requesting it raises `MissingRuleException`. A successful
planning run produces the ingest h5ad, `validation_report.json`,
`preprocessing_plan.json`, `qc_explore` artifacts, and `plan_review.awaiting_approval`.

S0 execution location is determined by the configured mode. **`run` is local-only;
`submit` is cluster-only** — Processing-MuAgent never submits cluster jobs itself.
Both refuse to start until execution mode is user-confirmed (Section 4) — so do not
reach this step before the user has chosen local vs HPC.

**HPC mode (`execution.mode` is `pbs` or `slurm`) — submit the planning head-job:**

```
source deliverables/plan/config/hpc.env
executor submit --config $CFG --executor pbs|slurm
executor hpc-status --config $CFG             # one-shot: report, then re-poll on a scheduled wakeup
```

Omit `--target` — `submit` auto-infers `plan_review_propose`. S0's QC exploration
needs 100+ GB, so it runs on a compute node inside the head-job (never the login
node). `submit` routes through Execution-MuAgent (kill-on-hang, survives SSH
disconnect) and returns in ≤90 s. After `submit`, follow **report-and-repoll** (see
`interaction_flow.md`): report one-shot `hpc-status`, then re-poll on a non-blocking
scheduled wakeup (~295s, per the `Next check:` line) — reporting only when the `State:`
fingerprint changes — until `monitor.pid` is removed or `plan_review` becomes
`awaiting_approval`.

**Local mode (`execution.mode` is `local`) — run everything locally:**

```
executor run --config $CFG --target plan_review_propose
```

Runs P1 → S0 (+ plan assembly + gate-arming). Small inputs: ~30s.

**If S0 OOMs:** in HPC mode, raise `PMA_RESOURCES_SCALE` (`configure-execution
--resources-scale N`) and `submit` again (no `--target`). In local mode, the machine
is too small — switch to HPC (`configure-execution --mode slurm|pbs`) and submit.
There is no automatic local→cluster retry.

Do **not** cluster-retry logic errors (pairing ambiguous, path missing, branch mismatch). Relay and let the user fix inputs or `declare-branch`.

## What to surface back

After `executor init`: confirm the canonical config path and blank context template path.

After biological-context write (cases a/b/c): confirm the file exists at `deliverables/plan/config/biological_context.md` and that you populated the fields the user told you about.

After `configure-execution`: confirm `execution.mode` **and `compute.device` (cpu/gpu)** and, for HPC, the path to `deliverables/plan/config/hpc.env`. Tell the user to `source` that file before any cluster submit/resume. When `device=gpu` on SLURM, also confirm the GPU partition/gres and the pinned `gpu_image_uri` (the container image is **pulled** from that registry reference — recorded in `site.config` / `~/.muagene/machine.config`, not a conda env).

After the planning phase completes (`plan_review_propose` — local `run` or HPC `submit`):

- If `deliverables/plan/context_summary.md` exists, paste its content back verbatim. Any conflicts (e.g. "report says mouse, file fingerprint says human") surface here and must be resolved before Step 3.
- If `executor run` errored:
  - **Phase 1 gate error** — context template is blank and user didn't opt out. Ask for context OR tell them to re-invoke with `--no-context`.
  - **S0 declared-vs-detected mismatch** — relay the raised error and ask the user to fix the declaration or the inputs.
  - **S0 ambiguous pairing** — relay the raised error; ask paired vs separate; re-run `executor declare-branch` and re-try.

Transition to Step 3 (`plan_review`) once `plan_review_<run>.md` exists under `deliverables/plan/` (written by `plan_review_propose` or `Processing-MuAgent plan-review`).

## Explicit non-actions

- Do NOT write `parameters.yaml` directly — use `executor declare-branch` / `executor revise`.
- Do NOT touch `biological_context.md` directly — always route through `context_mapper.write_report`.
- Do NOT copy input files into the run dir. Paths are referenced in place.
- Do NOT auto-retry on S0 errors. Relay the message; let the user correct the root cause.
- Do NOT skip biological context "to save time" — the Phase 1 gate exists because context shapes downstream QC thresholds; the cost of silently proceeding is recovering later from wrong defaults.
