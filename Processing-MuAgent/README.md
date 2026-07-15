# Processing-MuAgent

Processing-MuAgent guides reproducible preprocessing of single-cell RNA and ATAC
data. It combines a conversational scientific workflow with a deterministic pipeline,
pausing at the decisions that require biological judgment.

It produces quality-controlled, per-modality representations through clustering and
UMAP. It does **not** perform multimodal integration, cell-type annotation, marker
discovery, or gene-regulatory network inference.

Processing-MuAgent owns the scientific plan and user interaction. Its sibling
[Execution-MuAgent](../Execution-MuAgent/README.md) owns machine setup and cluster
execution.

## Supported analyses

| Workflow | Inputs | Result |
|---|---|---|
| `rna_only` | RNA counts | Processed RNA data |
| `atac_only` | ATAC fragments | Processed ATAC data |
| `paired` | RNA and ATAC from the same cells | A paired multimodal object with shared cells |
| `unpaired` | RNA and ATAC from different cell sets | Independent processed RNA and ATAC objects |

For paired data, barcode diagnostics verify that the modalities represent the same
cells. Direct overlap, subset relationships, common barcode suffixes, and an optional
translation table are considered. If pairing cannot be validated, the pipeline stops
and asks whether to correct the inputs, provide a translation, or explicitly continue
as unpaired data. It never changes the declared branch silently.

## Inputs

You provide:

- a run location for outputs;
- the intended workflow;
- a genome assembly;
- RNA and/or ATAC input locations;
- biological context such as organism, tissue, assay, and relevant publications;
- a confirmed execution mode: local or SLURM.

Supported RNA inputs:

| Format | Requirements |
|---|---|
| 10x HDF5 | Cell Ranger gene-expression or ARC matrix |
| 10x MEX | Matrix, barcode, and feature tables |
| AnnData | Raw integer counts in `.X` |
| Dense text | Gzipped, tab-delimited genes-by-cells matrix |

Supported ATAC inputs:

| Format | Requirements |
|---|---|
| Fragments | Block-gzipped fragments with a matching tabix index |
| BED4 | Gzipped four-column intervals; converted during intake |

Optional inputs include a raw RNA matrix for ambient-RNA estimation, user-defined
peaks, cell metadata, a barcode translation, and a biological-context document.

## How preprocessing works

```text
P1 context extraction
  → S0 ingest (load + validate + assemble plan + QC exploration and threshold preview)
  → (CHECKPOINT 1) plan_review
  → S1a ambient RNA correction → S1 RNA QC → S2 ATAC QC → S3 doublets
  → (CHECKPOINT 2) post_qc_review
  → qc_handoff (post-QC h5mu + integration manifest; canonical post-QC input for S4/S5)
  → user confirms the finish batch
  → S4 RNA normalization + HVG → S5 ATAC spectral embedding
  → S6 RNA PCA + per-modality neighbor graphs
  → S7 clustering (fixed resolutions: RNA 0.7, ATAC 0.5) → S8 UMAP
  → manifest
```

Stages that do not apply to an RNA-only or ATAC-only workflow are omitted automatically.

### 1. Extract context, ingest, and preview QC

P1 records the study context. S0 loads the data, checks input formats and genome
compatibility, verifies the declared workflow, and assembles the preprocessing plan.

Before plan review, S0 also explores the observed RNA and ATAC quality distributions.
It generates diagnostic plots, previews the effective QC thresholds, and estimates how
many cells each proposed filter would retain or remove. These are previews only; no QC
filter is applied until you approve the plan.

### 2. Review the preprocessing plan

Before QC runs, you review:

- the detected modalities and pairing assessment;
- ambient-RNA handling for RNA data;
- optional marker genes for comparing expression before and after ambient-RNA correction;
- proposed RNA and ATAC quality filters;
- doublet-detection strategy;
- dimensionality-reduction and clustering choices.

You can accept the proposal or revise supported parameters.
Plan-stage revisions are non-destructive because preprocessing has not started.

### 3. Run modality-aware QC

Recommended cutoffs are derived from the observed distributions where appropriate.
MAD-derived values are recalculated from each dataset during S0 exploration, and the
generated plan shows the effective cutoffs that will be applied. At either review
checkpoint, you can tune these policies, pin an exact MAD-derived bound, remove one
side of a bound, or disable an individual metric entirely.

**Default RNA cell filtering metrics**
For RNA, the pipeline can correct ambient RNA, filter cells using count, detected-gene,
mitochondrial, and ribosomal metrics, and detect doublets.

| Metric | Default policy |
|---|---|
| Total counts | Log-MAD bounds (`k = 5`), with a minimum lower floor of 500 |
| Detected genes | Log-MAD bounds (`k = 5`), with a minimum lower floor of 250 |
| Mitochondrial fraction | MAD-derived upper bound (`k = 3`), constrained to 5–20% |
| Ribosomal fraction | Maximum 50% |

**Default ATAC cell filtering metrics**
For ATAC, it evaluates fragment depth, transcription-start-site enrichment,
nucleosome signal, fraction of reads in peaks, and doublet probability.

| Metric | Default policy |
|---|---|
| Fragment count | Log-MAD bounds (`k = 5`), with a minimum lower floor of 1,000 |
| TSS enrichment | 1.5–50 |
| Nucleosome signal | Maximum 3 |
| Fraction of reads in peaks (FRiP) | Minimum 0.2 (20%) |

**Doublet filtering**

| Modality | Default score threshold |
|---|---|
| RNA (Scrublet) | 0.25 |
| ATAC (SnapATAC2) | 0.5 |

For paired data, cells flagged by either the RNA or ATAC detector are removed
(union policy), then only barcodes retained in both modalities continue.

### 4. Review observed QC

The second formal checkpoint summarizes cell retention, applied thresholds, doublet
calls, and diagnostic figures. You can:

- approve the observed QC;
- revise, pin, or completely disable individual QC metrics and rerun affected QC steps;
- complete or waive an optional marker-gene comparison for ambient-RNA correction;
- stop and resume later.

Post-QC revisions can invalidate existing QC results, so the agent previews their
impact and asks before applying them.

### 5. Create the handoff and finish

After QC approval, Processing-MuAgent creates a stable post-QC handoff for downstream
steps and verifies it. It then asks before starting the unattended finish batch.

Regenerable intermediates are cleaned as checkpoints complete, while durable QC and
handoff outputs are retained. Later QC revisions regenerate the affected steps.

The finish batch performs:

- RNA normalization, highly variable gene selection, PCA, and neighbor construction;
- ATAC feature selection, spectral representation, and neighbor construction;
- Leiden clustering and UMAP separately for each modality.

Paired workflows retain diagnostic RNA and ATAC clusterings; they do not create a joint
multimodal embedding. The agent stops after the processed data and manifests are ready.

## Your decision points

| When | Your action |
|---|---|
| Intake | Confirm workflow, inputs, context, and local or SLURM execution |
| Plan review | Accept or revise the proposed preprocessing strategy |
| QC review | Inspect observed QC; approve or revise |
| Post-QC handoff | Confirm whether to start the final unattended batch |

These are the only points that require scientific or operational input on the normal
path. The agent reports failures and asks for a decision instead of silently changing
data interpretation or execution settings.

## Outputs

The run produces four user-facing groups:

- **Plan:** biological context, validated inputs, chosen strategy, and provenance.
- **QC:** review summaries, figures, retained cells, and the post-QC integration handoff.
- **Results:** processed RNA and/or ATAC objects, per-modality embeddings, and a review
  notebook.
- **Manifests:** machine-readable records that describe the handoff and final outputs.

### Run directory layout

Comments identify each output; brackets show the generating pipeline stage.

```text
<run_dir>/
├── deliverables/
│   ├── plan/
│   │   ├── config/...                        # run configuration [intake]
│   │   ├── context_summary.md                # context summary [P1]
│   │   ├── plan_review_<run>.md              # plan-review document [plan_review]
│   │   └── plan_summary_<run>.html           # web plan review [plan_review]
│   ├── figures/                              # pipeline figures [S0–S8, as applicable]
│   ├── qc/
│   │   ├── qc_review_<run>.md                # QC-review document [post_qc_review]
│   │   ├── qc_summary_<run>.html             # web QC review [post_qc_review]
│   │   ├── post_qc_<run>.h5mu                # canonical post-QC handoff [qc_handoff]
│   │   ├── peaks_<run>.bed                   # per-sample peaks [qc_handoff; ATAC]
│   │   └── post_qc_manifest.json             # post-QC handoff manifest [qc_handoff]
│   └── results/
│       ├── processed_<run>.h5mu              # paired processed object [S8]
│       ├── rna_processed.h5ad                 # RNA processed object [S8]
│       ├── atac_processed.h5ad                # ATAC processed object [S8]
│       ├── review_processed_<run>.{ipynb,py}  # review notebook and script [manifest]
│       ├── run_manifest.json                  # preprocessing handoff manifest [manifest]
│       └── layout.json                        # deliverable layout record [manifest]
└── internal/...                               # pipeline-managed state [all stages]
```

Only files relevant to the selected workflow branch and execution mode are created.
Raw input files are referenced in place and never overwritten.

## Project context

Processing-MuAgent is MuAgene's scientific component. See the
[MuAgene guide](../README.md) for installation, framework-level responsibilities, and
how to start a run. The concise agent contract is documented in [AGENT.md](AGENT.md);
implementation procedures remain in the agent instructions rather than this user guide.
