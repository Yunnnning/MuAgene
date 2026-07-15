---
name: inputs_intake
domain: intake
purpose: Collect paths + biological context + execution mode, scaffold the run (init), and declare the branch. The SSOT for execution-mode intake heuristics.
activation: analysis type known; run.yaml not written yet
inputs: [user dialogue, raw input paths, genome assembly]
outputs: [deliverables/plan/config/run.yaml, deliverables/plan/config/biological_context.md, parameters.yaml, "SLURM only: hpc.env + site.config"]
calls_tools: [init, hpc-info, configure-execution, declare-branch]
reads_contracts: [run_yaml, site_config]
writes_state: [run.yaml, biological_context.md, parameters.yaml, "SLURM only: hpc.env + site.config"]
handoff: { next: plan_confirm, when: run scaffolded + branch declared + exec-mode confirmed, on_error: troubleshooting }
---

# Inputs intake — paths + optional biological context

Script for the turn(s) after the user has declared their analysis type
([`entry_declare.md`](entry_declare.md)). Goal: collect enough to build a valid `run.yaml`,
write it via `executor init`, populate biological context (if offered), configure execution
mode (local vs HPC), and declare the branch. Then hand off to
[`plan_confirm.md`](plan_confirm.md), which kicks off planning compute and drives the
plan_review gate.

## What to say

Tailor the required paths to the declared `workflow_branch`:

### For `rna_only`

> I need:
> - **RNA input path** — one of: 10x Cell Ranger `.h5`, 10x MEX directory (matrix.mtx + barcodes.tsv + features.tsv), or `.h5ad`.
> - **Genome assembly** — e.g. `mm10`, `GRCh38`. Required; no default. I'll cross-check it against features in the matrix where I can.
>
> Optional:
> - **Biological context** — organism, tissue, assay, any DOIs. Free text is fine, or paste a filled Biological Context Report, or give me a path to one.
> - **Seed** — default 42.
> - **Execution** — Should I run locally on this machine, or submit jobs to a cluster (HPC: SLURM)? - If you choose HPC and have not yet set the required `PMA_*` environment variables, I'll run `hpc-info` on the login node, list available queues/partitions, suggest a project code or account where I can detect one, and ask you to confirm before I write `hpc.env`.

### For `atac_only`

> I need:
> - **ATAC fragments path** — `fragments.tsv.gz`. The `.tbi` index must sit next to it with the same stem.
> - **Genome assembly** — e.g. `mm10`, `GRCh38`. Required. I cross-check the fragments' chromosome names against it (`chr1..chrY,chrM` vs `1..22,X,Y,MT`).
>
> Optional: same as rna_only.

### For `paired` or `unpaired`

> I need:
> - **RNA input path** (h5 / MEX / h5ad)
> - **ATAC fragments path** (`fragments.tsv.gz` + matching `.tbi`)
> - **Genome assembly**
>
> Optional: same as above.
>
> If the RNA and ATAC barcodes have Jaccard overlap ≥80%, or one modality's barcodes
> are ≥80% contained in the other (typical when cell counts differ before QC), pairing
> is validated. If you declared `paired` but validation fails, I will stop and ask
> before switching to `unpaired`. Ambiguous overlap also stops for resolution.

## Actions

Once the user answers, execute in order:

### 1. Draft `run.yaml`

Build an in-memory dict from the user's answers (omit fields not supplied):

```yaml
run_dir: <user's run_dir>
rna_path: <optional>
atac_fragments_path: <optional>
genome_assembly: <required>
# s1a_ambient_method: <auto | none | decontx | soupx>   # optional; else plan uses ingest
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
- `<run_dir>/deliverables/plan/` — created at init; `figures/`, `qc/`, and `results/` appear when outputs are written

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
- **Where to run** — `--mode local | slurm` (always ask; never auto-default).
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

**HPC (SLURM):**

1. Run `executor hpc-info`. Parse the JSON silently — do not dump the raw JSON to the user.
   Also read the **GPU** fields: `slurm.gpu_partitions`, `slurm.suggested_gpu_partition`,
   `slurm.suggested_gpu_gres` (SLURM). A non-empty `gpu_partitions` / a `suggested_gpu_gres`

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
   executor configure-execution --config $CFG --mode slurm \
       --slurm-partition <partition> --slurm-account <account> --confirmed-by-user
   ```
   `--confirmed-by-user` records the user's approval; without it `run`/`submit` refuse to launch.

   **If the user chose GPU**, add the device flags (sourced from `hpc-info`'s GPU fields):
   ```
   executor configure-execution --config $CFG --mode slurm \
       --slurm-partition <cpu_partition> --slurm-account <account> \
       --device gpu --gpu-partition <suggested_gpu_partition> --gpu-gres <suggested_gpu_gres> \
       --gpu-image-uri docker://<registry>/muagene-gpu:<tag> --confirmed-by-user
   ```
   `configure-execution` fails loud on missing prerequisites — pre-empt them: SLURM `--device gpu`
   **requires** `--gpu-gres` and `--gpu-image-uri` (the SLURM GPU env is a container PULLED from that
   pinned reference — or set `gpu_image_uri` once in `~/.muagene/machine.config` via init-machine);
   Add `--singularity-module <module>` when the site needs `module load` for singularity.

Do **not** invent partition/account names — use `hpc-info` results only. If `hpc-info` returns empty lists, ask the user for the values directly.

### 5. Declare the branch

```
executor declare-branch <rna_only|atac_only|paired|unpaired> --config $CFG
```

This writes `plan.workflow_branch_declared` to `parameters.yaml` as a `source=user`
assertion. S0 validates it against pairing diagnostics. A declared `paired` run that
cannot be validated stops with three choices: provide a barcode translation, correct
the inputs, or explicitly re-declare `unpaired`. No branch change is automatic.

## What to surface back

After `executor init`: confirm the canonical config path and blank context template path.

After biological-context write (cases a/b/c): confirm the file exists at `deliverables/plan/config/biological_context.md` and that you populated the fields the user told you about.

After `configure-execution`: confirm `execution.mode` **and `compute.device` (cpu/gpu)** and, for HPC, the path to `deliverables/plan/config/hpc.env`. Tell the user to `source` that file before any cluster submit/resume. When `device=gpu` on SLURM, also confirm the GPU partition/gres and the pinned `gpu_image_uri` (the container image is **pulled** from that registry reference — recorded in `site.config` / `~/.muagene/machine.config`, not a conda env).

Once the run is scaffolded, the branch is declared, and execution mode is confirmed, hand off
to [`plan_confirm.md`](plan_confirm.md) — it runs the planning phase (P1 → S0, which assembles the plan in-process), surfaces
`context_summary.md` / `validation_report.json`, and drives the plan_review gate. Any S0 or
context-gate errors at that point → [`troubleshooting.md`](troubleshooting.md).

## Explicit non-actions

- Do NOT write `parameters.yaml` directly — use `executor declare-branch` / `executor revise`.
- Do NOT touch `biological_context.md` directly — always route through `context_mapper.write_report`.
- Do NOT copy input files into the run dir. Paths are referenced in place.
- Do NOT auto-retry on S0 errors. Relay the message; let the user correct the root cause.
- Do NOT skip biological context "to save time" — the Phase 1 gate exists because context shapes downstream QC thresholds; the cost of silently proceeding is recovering later from wrong defaults.
