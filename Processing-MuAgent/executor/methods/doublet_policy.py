"""Doublet overlap categorisation + union-based removal policy."""
from __future__ import annotations

from typing import Any

import pandas as pd


def four_way_overlap(rna_flag: pd.Series, atac_flag: pd.Series) -> dict[str, int]:
    """rna_flag, atac_flag indexed by cell barcode (same index). Booleans."""
    rna_only = int(((rna_flag) & (~atac_flag)).sum())
    atac_only = int(((~rna_flag) & (atac_flag)).sum())
    both = int((rna_flag & atac_flag).sum())
    neither = int(((~rna_flag) & (~atac_flag)).sum())
    return {"rna_only": rna_only, "atac_only": atac_only, "both": both, "neither": neither}


def recommend_policy() -> dict[str, Any]:
    """Paired multiome always uses union: remove cells flagged by either detector.

    Doublet detectors are prone to false negatives; union is stricter and reduces
    contamination risk in downstream clustering.
    """
    return {
        "recommendation": "union",
        "rationale": (
            "Paired multiome uses union policy: remove cells flagged by EITHER detector "
            "(RNA Scrublet or ATAC SnapATAC2). Detectors are prone to false negatives, "
            "so union minimises doublet contamination."
        ),
    }


def apply_policy(rna_flag: pd.Series, atac_flag: pd.Series, policy: str) -> pd.Series:
    if policy != "union":
        raise ValueError(f"policy must be 'union', got {policy!r}")
    return rna_flag | atac_flag
