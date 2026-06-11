"""Shared per-cell RNA QC-metric computation.

Single implementation of the mitochondrial / ribosomal gene flagging and
``scanpy`` per-cell QC metric computation used by both the pre-plan QC
exploration (``executor.qc_explore``) and the real RNA QC stage
(``executor.stages.s1_rna_qc``). The two callers run on different inputs
(pre- vs post-ambient RNA) so the *compute* is intentionally not shared, but
the metric definition is, so the preview and the stage agree on what each
metric means.
"""
from __future__ import annotations

from typing import Any

# Symbols vary by species/case (Mt-/mt-/MT-/mt: for fly, etc.) and by whether
# var_names are gene IDs (Ensembl) or symbols. Use a permissive case-insensitive
# prefix set.
MT_PATTERN = r"(?i)^(?:mt[-:_]|mito[-:_]?)"      # mt-, MT-, Mt:, mito-, etc.
RIBO_PATTERN = r"(?i)^(?:rps|rpl|mrps|mrpl)"     # cytoplasmic + mito ribosomal proteins


def compute_rna_qc_metrics(adata: Any) -> dict[str, int]:
    """Flag mt/ribo genes and compute scanpy per-cell QC metrics in place.

    Adds ``adata.var["mt"]`` / ``adata.var["ribo"]`` and the standard
    ``calculate_qc_metrics`` columns (``total_counts``, ``n_genes_by_counts``,
    ``pct_counts_mt``, ``pct_counts_ribo``). Returns the detected mt/ribo gene
    counts so callers can sanity-check / log.
    """
    import scanpy as sc

    adata.var["mt"] = adata.var_names.astype(str).str.contains(MT_PATTERN, regex=True, na=False)
    adata.var["ribo"] = adata.var_names.astype(str).str.contains(RIBO_PATTERN, regex=True, na=False)
    n_mt = int(adata.var["mt"].sum())
    n_ribo = int(adata.var["ribo"].sum())
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None,
                               log1p=False, inplace=True)
    return {"n_mt": n_mt, "n_ribo": n_ribo}
