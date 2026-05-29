"""Doublet overlap categorisation + goal-based recommendation."""
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


def recommend_policy(study_goal: str | None) -> dict[str, Any]:
    """study_goal: 'clustering_inference' | 'rare_populations' | None.

    Default (clustering_inference / unspecified) → union (remove cells flagged
    by EITHER detector). intersection only when study_goal=rare_populations.
    """
    g = (study_goal or "").strip().lower()
    if g == "rare_populations":
        return {
            "recommendation": "intersection",
            "rationale": "study_goal=rare_populations prioritises avoiding false removals; remove only cells flagged by BOTH detectors.",
        }
    if g == "clustering_inference" or not g:
        return {
            "recommendation": "union",
            "rationale": (
                "study_goal=clustering_inference (or unspecified fallback) prioritises minimizing "
                "doublet contamination; remove cells flagged by EITHER detector."
            ),
        }
    return {
        "recommendation": "union",
        "rationale": f"Unknown study_goal={study_goal!r}; default fallback is UNION.",
    }


def apply_policy(rna_flag: pd.Series, atac_flag: pd.Series, policy: str) -> pd.Series:
    if policy == "union":
        return rna_flag | atac_flag
    if policy == "intersection":
        return rna_flag & atac_flag
    raise ValueError(f"policy must be 'union' or 'intersection', got {policy!r}")
