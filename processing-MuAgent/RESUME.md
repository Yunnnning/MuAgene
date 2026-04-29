# Resume checkpoint — 2026-04-23

## Status

Implementation of `processing-MuAgent` is complete. End-to-end run on the example
data has been **deliberately deferred** at the user's request.

## What's been built

- Directory structure scaffolded under `/Users/yunning/Desktop/MuGene-agents/processing-MuAgent/`.
- **Executor modules** in `executor/`:
  - `hashing.py`, `provenance.py`, `log.py`, `approval.py`
  - `io.py` (RNA format autodetect + fragment validation + genome fingerprint)
  - `pairing.py` (four-strategy paired/separate/ambiguous detection)
  - `metadata.py` (minimal reconstruction + recovery)
  - `context.py` (Biological Context Report parser + deterministic inference)
  - `doi_fetch.py` (Crossref lookup + per-DOI cache)
  - `plan_assembler.py` (P2 `preprocessing_plan.json`)
  - `manifest.py` (`run_manifest.json`, handoff contract v1.0.0)
- **Method helpers** in `executor/methods/`:
  - `mad_thresholds.py` (log-MAD bounds, upper-bound helper)
  - `resolution_sweep.py` (Leiden sweep + silhouette + seed stability + stable-knee)
  - `doublet_policy.py` (four-way overlap + goal-based recommendation)
- **Stage scripts** in `executor/stages/` — S0 ingest, S1 RNA QC, S2 ATAC QC,
  S3 doublets, S4 RNA norm + HVG, S5 ATAC TF-IDF + LSI, S6 dim reduction + neighbors,
  S7 clustering (sweep-based resolution), S8 UMAP (paired h5mu / separate h5ad).
- **Snakemake workflow** in `workflow/`:
  - `Snakefile` + per-stage `rules/*.smk` implementing `<stage>_propose` + `<stage>_execute`
    rule pairs, each execute rule gated on `checkpoints/<stage>.approved`.
  - `rules/manifest.smk` assembles the final `run_manifest.json`.
- **CLI** (`executor/cli.py`): `init | propose | approve | revise | status | run`.
- **Config example**: `config/run.example.yaml` pointed at the example data.
- **README.md** with usage + scope.

## What's NOT done yet

- **End-to-end run on the example data.** Not attempted yet per user pause.
- **Debug pass** to fix whatever breaks on first run. Expected hotspots listed below.
- **Unit tests** scaffolded in the plan but not yet written.
- **Per-rule conda env YAMLs** under `workflow/envs/` (the implementation runs
  against the shared `cell_annotation` env for MVP).

## Environment (already set up)

- Python: `/Users/yunning/micromamba/envs/cell_annotation/bin/python` (3.10.18)
- Pre-installed: scanpy 1.11.4, anndata 0.10.9, snapatac2 2.8.0, umap-learn
- Added (via pip into `cell_annotation`): muon 0.1.7, mudata 0.3.6, scrublet,
  leidenalg, snakemake 7.32.4

## Exact commands to resume next session

```bash
PY=/Users/yunning/micromamba/envs/cell_annotation/bin/python
PROJ=/Users/yunning/Desktop/MuGene-agents/processing-MuAgent

# 1. Install the executor package in editable mode (first-time only)
$PY -m pip install -e "$PROJ"

# 2. Copy the example config and initialise the run directory
cp "$PROJ/config/run.example.yaml" /Users/yunning/Desktop/MuGene-agents/example/output/run_example.yaml
$PY -m executor.cli init --config /Users/yunning/Desktop/MuGene-agents/example/output/run_example.yaml

# 3. (Optional) fill in the Biological Context Report
#   edit /Users/yunning/Desktop/MuGene-agents/example/output/run_example/biological_context.md
#   Suggested contents for this dataset (GSE268104, mouse brain multiome):
#     - Organism: Mus musculus (C57BL/6)
#     - Tissue / sample: (see GSE268104 description)
#     - Assay: 10x Genomics multiome snRNA-seq + snATAC-seq
#     - DOI(s): (optional; any related GSE268104 paper)

# 4. Run the pipeline end-to-end with auto-approval (noninteractive MVP path)
$PY -m executor.cli run \
  --config /Users/yunning/Desktop/MuGene-agents/example/output/run_example/run.yaml \
  --auto-approve

# 5. Inspect outputs
ls /Users/yunning/Desktop/MuGene-agents/example/output/run_example/artifacts/s8_umap/
cat /Users/yunning/Desktop/MuGene-agents/example/output/run_example/run_manifest.json
```

## Expected hotspots / likely first-run failures to debug

Listed in order of decreasing probability. These reflect MVP shortcuts — not
architecture changes. Each is a normal debug-loop iteration, not a redesign.

1. **SnapATAC2 API drift** — `s2_atac_qc.py` and `s5_atac_lsi.py` call
   `snap.pp.import_fragments`, `snap.metrics.tsse`, `snap.metrics.frag_size_distr`,
   `snap.pp.add_tile_matrix`, `snap.pp.select_features`, `snap.tl.spectral`,
   `snap.tl.leiden`, `snap.tl.umap`. Names vary by SnapATAC2 version (2.6 vs 2.7 vs
   2.8); the installed version is 2.8.0. First failure will likely be an
   `AttributeError` on one of these; look up the installed version's API and
   adapt. Genome helper `snap.genome.mm10` likewise may live under a different
   path in 2.8.
2. **SnapATAC2 AnnData subsetting** — `s0/s2/s3/s5` call `adata.subset(obs_indices=..., out=...)`
   and `adata.close()`. The exact API may differ; may need `snap.read`/`snap.write`
   patterns instead of the in-place object handle semantics.
3. **Scrublet on small / filtered data** — `s3_doublets.py` calls Scrublet;
   sometimes it fails to find a threshold on small datasets. A fallback
   (`scores > 0.2`) is already in place; may need to adjust if the example has
   very few cells surviving QC.
4. **h5 feature_types column name** — `io.detect_peaks_in_10x_h5` reads
   `/matrix/features/feature_type` directly from h5py. Different Cell Ranger
   versions name this differently (`feature_type` vs `_all_tag_keys`). Verify
   against the actual example h5.
5. **AnnData handle around SnapATAC2 writes** — `s5_atac_lsi.py` calls `adata.close()`
   before writing a separate summary. The s5 output h5ad may need to be produced
   via a different mechanism (`adata.to_memory()` then `write_h5ad`).
6. **UMAP single-threaded** — `NUMBA_NUM_THREADS` / `OMP_NUM_THREADS` not yet set
   in the CLI `run` path. If results are non-deterministic across re-runs, wrap
   the snakemake subprocess env with these.

## Plan-file reference

The approved design is at `/Users/yunning/.claude/plans/you-are-tasked-with-zippy-moonbeam.md`.
If any pipeline behaviour diverges from the plan, treat the plan as authoritative
and adjust the implementation.

## Task list at pause

- #1 Scaffold — ✅ completed
- #2 Env + deps — ✅ completed
- #3 Executor modules — ✅ completed
- #4 P1 + P2 — ✅ completed
- #5 Stages S0–S8 — ✅ completed
- #6 Snakefile + CLI — ✅ completed
- **#7 Run end-to-end + debug — PAUSED** (resume here)
