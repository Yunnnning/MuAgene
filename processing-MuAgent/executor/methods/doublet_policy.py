"""Doublet overlap categorisation + branch-aware removal.

Removal rule (no longer configurable):
    - paired       : intersection — remove cells flagged by BOTH detectors.
    - rna_only     : remove cells flagged by Scrublet (RNA).
    - atac_only    : remove cells flagged by the ATAC scrublet detector.
    - separate     : remove per-cell whatever detector flagged it (the two
                     modalities live on disjoint barcode sets here, so
                     intersection would be empty).

The previous `union` policy and `study_goal`-driven recommendation have been
removed. Raw per-cell scores + boolean flags from each detector are still
preserved in `calls.parquet` for reproducibility.
"""
from __future__ import annotations

import pandas as pd


def four_way_overlap(rna_flag: pd.Series, atac_flag: pd.Series) -> dict[str, int]:
    """rna_flag, atac_flag indexed by cell barcode (same index). Booleans."""
    rna_only = int(((rna_flag) & (~atac_flag)).sum())
    atac_only = int(((~rna_flag) & (atac_flag)).sum())
    both = int((rna_flag & atac_flag).sum())
    neither = int(((~rna_flag) & (~atac_flag)).sum())
    return {"rna_only": rna_only, "atac_only": atac_only, "both": both, "neither": neither}


def combine_flags(rna_flag: pd.Series, atac_flag: pd.Series,
                  workflow_branch: str) -> pd.Series:
    """Return the boolean removal mask for the merged per-barcode table.

    Inputs are already NaN-filled (False for missing modality) by the caller.
    """
    if workflow_branch == "paired":
        return rna_flag & atac_flag
    # rna_only / atac_only / separate: cells live on disjoint barcode sets
    # across modalities, so the right operation is "remove what was flagged by
    # whichever detector saw this barcode" — equivalent to OR with all-False
    # placeholders for the missing modality.
    return rna_flag | atac_flag


def removal_rule(workflow_branch: str) -> str:
    """Human-readable description of the removal rule for a given branch.
    Used in provenance and QC summary."""
    if workflow_branch == "paired":
        return ("intersection — remove only cells flagged by BOTH Scrublet (RNA) "
                "and SnapATAC2 scrublet (ATAC); preserves modality-specific signal.")
    if workflow_branch == "rna_only":
        return "Scrublet only (RNA detector); single-modality run."
    if workflow_branch == "atac_only":
        return "SnapATAC2 scrublet only (ATAC detector); single-modality run."
    return ("per-cell available detector (separate branch: RNA and ATAC live on "
            "disjoint barcode sets so intersection is not meaningful).")
