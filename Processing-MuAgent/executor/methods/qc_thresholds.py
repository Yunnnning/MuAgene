"""Shared QC threshold derivation.

Single source of truth for the MAD-based threshold math used by both the real QC
stages (``s1_rna_qc`` / ``s2_atac_qc``) and the pre-plan QC exploration
(``executor.qc_explore``). Keeping the derivation here guarantees the exploration
preview matches exactly what the stages will compute. Pure math only — no
provenance writes, no I/O.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from . import mad_thresholds as _mad


# --- RNA (S1) --------------------------------------------------------------

def rna_thresholds(
    obs: Any,
    *,
    total_counts_k_mad: float,
    n_genes_k_mad: float,
    pct_mt_k: float,
    pct_mt_ceiling: float,
    pct_mt_floor: float,
    min_counts_floor: float,
    min_genes_floor: float,
) -> dict[str, float]:
    """Derive RNA QC bounds from a QC-metrics ``obs`` frame (mirrors S1).

    Requires columns ``total_counts``, ``n_genes_by_counts``, ``pct_counts_mt``.
    MAD bounds are computed on the count-floor-passing subset, then lower bounds
    are clamped to the absolute floors.
    """
    tc = np.asarray(obs["total_counts"], dtype=float)
    ng = np.asarray(obs["n_genes_by_counts"], dtype=float)
    mt = np.asarray(obs["pct_counts_mt"], dtype=float)

    keep_floor = tc >= float(min_counts_floor)
    c_lo, c_hi = _mad.log_mad_bounds(tc[keep_floor], k=total_counts_k_mad)
    c_lo_mad_raw = float(c_lo)
    c_lo = max(c_lo, float(min_counts_floor))
    g_lo, g_hi = _mad.log_mad_bounds(ng[keep_floor], k=n_genes_k_mad)
    g_lo_mad_raw = float(g_lo)
    g_lo = max(g_lo, float(min_genes_floor))
    mt_subset = mt[keep_floor]
    pct_mt_mad_raw = _mad.mad_upper_raw(mt_subset, k=pct_mt_k)
    pct_mt_upper = _mad.upper_bound(
        mt_subset, k=pct_mt_k, floor=pct_mt_floor, ceiling=pct_mt_ceiling
    )
    return {
        "total_counts_min": float(c_lo),
        "total_counts_max": float(c_hi),
        "total_counts_mad_lo_raw": c_lo_mad_raw,
        "n_genes_min": float(g_lo),
        "n_genes_max": float(g_hi),
        "n_genes_mad_lo_raw": g_lo_mad_raw,
        "pct_counts_mt_max": float(pct_mt_upper),
        "pct_counts_mt_mad_raw": float(pct_mt_mad_raw),
    }


def rna_pass_masks(
    obs: Any, th: dict[str, float], *, pct_ribo_max: float
) -> dict[str, np.ndarray]:
    """Per-metric boolean *pass* masks on the full ``obs`` (mirrors S1)."""
    tc = np.asarray(obs["total_counts"], dtype=float)
    ng = np.asarray(obs["n_genes_by_counts"], dtype=float)
    mt = np.asarray(obs["pct_counts_mt"], dtype=float)
    ribo = np.asarray(obs["pct_counts_ribo"], dtype=float)
    return {
        "total_counts": (tc >= th["total_counts_min"]) & (tc <= th["total_counts_max"]),
        "n_genes": (ng >= th["n_genes_min"]) & (ng <= th["n_genes_max"]),
        "pct_counts_mt": mt <= th["pct_counts_mt_max"],
        "pct_counts_ribo": ribo <= float(pct_ribo_max),
    }


# --- ATAC (S2) -------------------------------------------------------------

def atac_n_fragment_bounds(
    n_frag: np.ndarray, *, k_mad: float, n_frag_floor: float
) -> tuple[float, float, float]:
    """MAD bounds on log(n_fragments) after the absolute floor (mirrors S2).

    Returns ``(applied_lower, upper, mad_lower_raw)`` where ``mad_lower_raw`` is
    the log-MAD lower bound before the absolute ``n_frag_floor`` clamp.
    """
    n_frag = np.asarray(n_frag, dtype=float)
    f_lo_mad_raw = float(n_frag_floor)
    if n_frag.size:
        keep_floor = n_frag >= float(n_frag_floor)
        if keep_floor.any():
            f_lo, f_hi = _mad.log_mad_bounds(n_frag[keep_floor], k=k_mad)
            f_lo_mad_raw = float(f_lo)
        else:
            f_lo, f_hi = float(n_frag_floor), float(n_frag.max() if n_frag.size else 1e6)
    else:
        f_lo, f_hi = float(n_frag_floor), 1e12
    f_lo = max(f_lo, float(n_frag_floor))
    return float(f_lo), float(f_hi), f_lo_mad_raw


def atac_pass_masks(
    n_frag: np.ndarray,
    tss: np.ndarray,
    ns: np.ndarray,
    *,
    f_lo: float,
    f_hi: float,
    tss_min: float,
    tss_max: float,
    nuc_max: float,
    n_pre: int,
) -> dict[str, np.ndarray]:
    """Per-metric boolean *pass* masks for ATAC (mirrors S2; FRiP excluded)."""
    n_frag = np.asarray(n_frag, dtype=float)
    tss = np.asarray(tss, dtype=float)
    ns = np.asarray(ns, dtype=float)
    pass_frag = (
        (n_frag >= f_lo) & (n_frag <= f_hi) if n_frag.size
        else np.ones(n_pre, dtype=bool)
    )
    pass_tss = (
        (tss > tss_min) & (tss < tss_max) if tss.size
        else np.ones(n_pre, dtype=bool)
    )
    pass_ns = (
        ns < nuc_max if ns.size and np.isfinite(ns).any()
        else np.ones(n_pre, dtype=bool)
    )
    return {
        "n_fragments": pass_frag,
        "tss_enrichment": pass_tss,
        "nucleosome_signal": pass_ns,
    }
