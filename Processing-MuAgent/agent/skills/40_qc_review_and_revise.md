---
name: qc_review_and_revise
domain: QC
purpose: Drive the post_qc_review gate (#2), safely revise QC thresholds, close the marker-gene decision, and build the post-QC handoff.
activation: status shows post_qc_review awaiting_approval; remains active through qc_handoff verification and finish-batch confirmation
inputs: [deliverables/qc/qc_review_<run>.md, internal/proposals/post_qc_review.yaml, internal/parameters.yaml]
outputs: [internal/parameters.yaml, post_qc_review.approved, deliverables/qc/post_qc_<run>.h5mu, deliverables/qc/peaks_<run>.bed, deliverables/qc/post_qc_manifest.json]
calls_tools: [status, revise, "revise --dry-run", approve, marker-gene-check, propose, submit, run]
reads_contracts: [parameters, latest_snapshot]
writes_state: [parameters.yaml, post_qc_review.approved]
handoff: { next: run_execution, when: qc_handoff verified + user confirmed finish batch, on_error: troubleshooting }
---

# QC review and revise — post_qc_review gate (#2)

Relay `deliverables/qc/qc_review_<run>.md` **verbatim** and point the user to
`deliverables/qc/qc_summary_<run>.html`. This skill owns post-QC revision, approval,
`qc_handoff`, and confirmation before the finish batch. Plan-review revision is
non-destructive and stays in [`20_plan_confirm.md`](20_plan_confirm.md).

## Revise guardrail (post-QC only)

At this gate, `executor revise` deletes stale stage, downstream, gate, and prior handoff
outputs so Snakemake recomputes them. Before a real revision:

1. Diagnose the binding constraint. A lower bound is `max(MAD, floor)` and an upper
   bound is `min(MAD, ceiling)`; changing a non-binding knob has no effect. When median
   mitochondrial percentage is low, `pct_mt_floor` commonly binds instead of
   `pct_mt_ceiling`.
2. Preview the exact mutation and deletion set:
   ```bash
   executor revise <stage> <key>=<value> --config $CFG --dry-run
   ```
3. Confirm with the user, then apply the same command without `--dry-run`, adding
   `--rationale "<user's reason>"`.

Revise input knobs, not computed MAD output keys. Pin an exact MAD-derived cutoff with
its `*_override` key. Never edit state files or checkpoint sentinels, and do not hand-delete
artifacts; `revise` owns invalidation. If cleanup already removed QC h5ads,
the executor also clears the required S1/S1a markers on RNA branches.

If a revision was accidental, rerun the invalidated stages; there is no in-place undo.
See [`90_troubleshooting.md`](90_troubleshooting.md) for recovery.

## Re-run after revision

Revisions to `s1_rna_qc`, `s2_atac_qc`, or `s3_doublets` delete their required durable
outputs, so the automated stages rerun without per-stage approval commands. Do not use `--auto-approve`;
it can prematurely recreate the QC gate approval and trigger cleanup.

### SLURM

1. Check for an older head job for this run and cancel it if still active:
   ```bash
   squeue -u $USER | grep "pma_head_job_$(basename <run_dir>)"
   # scancel <JOBID>
   ```
2. Submit the invalidated QC path:
   ```bash
   source deliverables/plan/config/hpc.env
   executor submit --config $CFG --executor slurm
   executor hpc-status --config $CFG
   ```
3. Follow [`80_hpc_monitoring.md`](80_hpc_monitoring.md). The inferred
   `post_qc_review_propose` target reruns QC and regenerates the reports. Use
   `executor propose post_qc_review --config $CFG` only for a manual report refresh.

### Local

```bash
executor run --config $CFG --target s3_doublets_execute
executor propose post_qc_review --config $CFG
```

After either path, read `qc_review_<run>.md` and relay it **verbatim**. Do not substitute
`internal/proposals/post_qc_review.yaml`. Ask the user to approve QC or revise again.

## Marker-gene decision

If the QC report says **"Marker gene expression check not performed"**, relay that notice
and obtain an explicit decision before approval: provide genes and run the check, or
explicitly decline.

**Hard rule — never pick genes:** do not select, suggest, look up, or test gene names.
Ask the user to provide the symbols, then run:

```bash
executor marker-gene-check --config $CFG <gene1> <gene2> ...
```

The command refreshes the QC reports; `--plot-only` is for layout iteration. If
`internal/artifacts/s1a_ambient/tsne_coords_cache.parquet` exists and the cell set is
unchanged, run on the login node. Without the cache, use a memory-appropriate cluster job
for HPC runs or run inline locally. Complete or explicitly waive this check before QC
approval because approval cleans the source working h5ad.

## Approve QC and build the handoff

When thresholds are accepted and the marker-gene decision is closed:

1. Approve the human gate:
   ```bash
   executor approve post_qc_review --config $CFG
   ```
2. Immediately build `qc_handoff`:
   ```bash
   # SLURM
   source deliverables/plan/config/hpc.env
   executor submit --config $CFG --executor slurm --target qc_handoff
   executor hpc-status --config $CFG
   ```
   ```bash
   # Local
   executor run --config $CFG --target qc_handoff
   ```
3. Verify:
   - `deliverables/qc/post_qc_<run>.h5mu`
   - `deliverables/qc/peaks_<run>.bed` for ATAC branches
   - `deliverables/qc/post_qc_manifest.json`
4. Tell the user the handoff is ready and obtain explicit user approval before starting
   S4–S8. Do not submit the finish batch yet.

On paired runs, confirm the union doublet-removal policy here; it is not a separate gate.

## QC revision reference

### Honored input keys

| Stage | Keys |
|---|---|
| `s1_rna_qc` | `total_counts_k_mad`, `n_genes_k_mad`, `pct_mt_k`, `pct_mt_ceiling`, `pct_mt_floor`, `pct_ribo_max`, `min_counts_floor`, `min_genes_floor`, `min_cells_per_gene` |
| `s2_atac_qc` | `frip_min`, `tss_enrichment_min`, `tss_enrichment_max`, `nucleosome_signal_max`, `n_fragments_k_mad`, `n_fragments_floor` |
| `s3_doublets` | `rna_doublet_score_threshold`, `atac_doublet_probability_threshold` |

### Pin an exact MAD-derived bound

| Stage | Override keys |
|---|---|
| `s1_rna_qc` | `total_counts_min_override`, `total_counts_max_override`, `n_genes_min_override`, `n_genes_max_override`, `pct_counts_mt_max_override` |
| `s2_atac_qc` | `n_fragments_min_override`, `n_fragments_max_override` |

Example:

```bash
executor revise s1_rna_qc n_genes_min_override=300 --config $CFG --rationale "user requested at least 300 genes"
```

An override is applied even when it is more permissive than the recommended
floor/ceiling; the report and event log make that deviation visible. Fixed thresholds
such as `pct_ribo_max`, `frip_min`, TSS/nucleosome limits, floors/ceilings, and doublet
thresholds do not use an `_override` suffix.

### Skip individual metrics

RNA:

| User intent | Keys to set |
|---|---|
| Skip `total_counts` | `total_counts_k_mad=999`, `min_counts_floor=0` |
| Remove only the `total_counts` upper bound | `total_counts_k_mad=999` |
| Skip `n_genes` | `n_genes_k_mad=999`, `min_genes_floor=0` |
| Remove only the `n_genes` upper bound | `n_genes_k_mad=999` |
| Skip `pct_counts_mt` | `pct_mt_k=999`, `pct_mt_ceiling=100` |
| Skip `pct_counts_ribo` | `pct_ribo_max=100` |

ATAC:

| User intent | Keys to set |
|---|---|
| Skip `n_fragments` | `n_fragments_k_mad=999`, `n_fragments_floor=0` |
| Remove only the `n_fragments` upper bound | `n_fragments_k_mad=999` |
| Skip `tss_enrichment` | `tss_enrichment_min=0`, `tss_enrichment_max=999` |
| Remove only the `tss_enrichment` upper bound | `tss_enrichment_max=999` |
| Skip `nucleosome_signal` | `nucleosome_signal_max=999` |
| Skip `frip` | `frip_min=0` |

`*_k_mad` is symmetric: it moves both lower and upper MAD bounds. To restore a finite
upper MAD bound while keeping an exact lower floor, also set the corresponding
`*_min_override` to that floor.
