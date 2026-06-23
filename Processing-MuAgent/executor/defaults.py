"""Single source of truth for QC default parameter values.

Every QC threshold default lives here exactly once. The places that previously
each carried their own copy now read from this module, so they can never silently
drift:

- ``plan_assembler.assemble_plan`` — writes the defaults into
  ``preprocessing_plan.json`` (the authoritative plan layer).
- ``executor.figures`` — the ``DEFAULT_*`` reference constants used to draw the
  grey "default" markers; re-exported from here.
- ``executor.qc_explore`` — the pre-plan QC preview (via ``figures``' constants),
  so the preview matches what the stages compute.
- the per-stage ``_resolve_param(..., key, <default>)`` fallbacks in
  ``stages/s1_rna_qc``, ``stages/s2_atac_qc``, ``stages/s3_doublets``,
  ``stages/s7_clustering``.

Only *values* live here; the rationale / source / confidence prose stays in
``plan_assembler`` because it is plan-presentation text, not a tunable constant.
Types are the exact types the plan serialises (int floors, float multipliers) so
``preprocessing_plan.json`` stays byte-identical to the pre-centralisation output.

Pure data — no imports from the rest of ``executor`` (this is a leaf module so
the lightweight ``plan_assembler`` / stage modules can import it without pulling
in matplotlib via ``figures``).
"""
from __future__ import annotations

from typing import Any

# stage id -> {parameter name: default value}. Keep the keys/types in lockstep
# with plan_assembler's `p(...)` calls — test_harness_consistency enforces it.
QC_DEFAULTS: dict[str, dict[str, Any]] = {
    "s1_rna_qc": {
        "total_counts_k_mad": 5.0,
        "n_genes_k_mad": 5.0,
        "pct_mt_k": 3.0,
        "pct_mt_ceiling": 20.0,
        "pct_mt_floor": 5.0,
        "pct_ribo_max": 50.0,
        "min_cells_per_gene": 3,
        "min_counts_floor": 500,
        "min_genes_floor": 250,
    },
    "s2_atac_qc": {
        "tss_enrichment_min": 1.5,
        "tss_enrichment_max": 50.0,
        "n_fragments_k_mad": 5.0,
        "n_fragments_floor": 1000,
        "nucleosome_signal_max": 3.0,
        "frip_min": 0.2,
    },
    "s3_doublets": {
        "rna_doublet_score_threshold": 0.25,
        "atac_doublet_probability_threshold": 0.5,
    },
    "s7_clustering": {
        "rna_resolution": 0.7,
        "atac_resolution": 0.5,
        "random_state": 0,
    },
}

# Named reference constants consumed by figures.py (and, through it, qc_explore.py)
# for drawing the grey "default" markers. Single home here; figures.py re-exports.
# Floors are floats here (marker geometry) while the plan/stage copies stay int —
# the values are equal, only the rendering uses the float form.
DEFAULT_TOTAL_COUNTS_K_MAD = QC_DEFAULTS["s1_rna_qc"]["total_counts_k_mad"]
DEFAULT_N_GENES_K_MAD = QC_DEFAULTS["s1_rna_qc"]["n_genes_k_mad"]
DEFAULT_PCT_MT_K = QC_DEFAULTS["s1_rna_qc"]["pct_mt_k"]
DEFAULT_PCT_MT_CEILING = QC_DEFAULTS["s1_rna_qc"]["pct_mt_ceiling"]
DEFAULT_PCT_MT_FLOOR = QC_DEFAULTS["s1_rna_qc"]["pct_mt_floor"]
DEFAULT_PCT_RIBO_MAX = QC_DEFAULTS["s1_rna_qc"]["pct_ribo_max"]
DEFAULT_MIN_COUNTS_FLOOR = float(QC_DEFAULTS["s1_rna_qc"]["min_counts_floor"])
DEFAULT_MIN_GENES_FLOOR = float(QC_DEFAULTS["s1_rna_qc"]["min_genes_floor"])
DEFAULT_N_FRAG_K_MAD = QC_DEFAULTS["s2_atac_qc"]["n_fragments_k_mad"]
DEFAULT_N_FRAG_FLOOR = float(QC_DEFAULTS["s2_atac_qc"]["n_fragments_floor"])
DEFAULT_TSS_MIN = QC_DEFAULTS["s2_atac_qc"]["tss_enrichment_min"]
DEFAULT_TSS_MAX = QC_DEFAULTS["s2_atac_qc"]["tss_enrichment_max"]
DEFAULT_NUC_MAX = QC_DEFAULTS["s2_atac_qc"]["nucleosome_signal_max"]
