---
name: completion_handoff
domain: finish
purpose: After manifest, validate outputs, surface results + the Integration handoff bundle, then HARD STOP at S8.
activation: manifest stage complete (finish batch done)
inputs: [deliverables/results/run_manifest.json, deliverables/qc/post_qc_manifest.json]
outputs: []
calls_tools: [status, finish-cleanup]
reads_contracts: [run_manifest, post_qc_manifest]
writes_state: []
handoff: { next: STOP, when: outputs surfaced, on_error: troubleshooting }
---

# Completion + handoff — STOP at S8

Entered when `manifest` completes. Preprocessing is done; your job is to validate, surface,
and stop.

## Steps

1. Read `deliverables/results/run_manifest.json`; extract `workflow_branch` and `outputs`.
   Confirm each listed output path exists.
2. **Reclaim disk — `executor finish-cleanup --config <run.yaml>`.** Only after step 1
   confirms the processed deliverable exists and is valid. The command re-validates it,
   removes only regenerable S4–S8 working files, and preserves durable markers and
   deliverables. It refuses cleanup if validation fails.
3. Point the user at:
   - the processed data + `run_manifest.json` under `deliverables/results/`,
   - the review notebook `review_processed_<run>.ipynb` (load + inspect + re-cluster at a
     custom resolution),
   - the UMAP figures in `deliverables/figures/`,
   - the QC summary `deliverables/qc/qc_review_<run>.md`,
   - the **Integration handoff bundle** under `deliverables/qc/`:
     `deliverables/qc/post_qc_manifest.json` (`muagene.post_qc_handoff/1`) +
     `deliverables/qc/post_qc_<run>.h5mu` +
     `deliverables/qc/peaks_<run>.bed` (ATAC branches with a peak set).
4. One-line sign-off, then **stop**:
   > Run complete. Outputs at `deliverables/results/`. I stop here — integration, annotation,
   > marker discovery, and GRN are out of scope (different subagents).
