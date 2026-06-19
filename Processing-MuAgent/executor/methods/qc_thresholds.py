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


def _eff(derived: float, override: float | None) -> float:
    """Effective bound: a user override wins over the MAD/floor-derived value.

    The derived value is always computed and returned separately (the ``*_derived``
    keys) so the figures can still draw the MAD/fixed line in grey while the chosen
    override is drawn in red.
    """
    return float(override) if override is not None else float(derived)


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
    total_counts_min_override: float | None = None,
    total_counts_max_override: float | None = None,
    n_genes_min_override: float | None = None,
    n_genes_max_override: float | None = None,
    pct_counts_mt_max_override: float | None = None,
) -> dict[str, float]:
    """Derive RNA QC bounds from a QC-metrics ``obs`` frame (mirrors S1).

    Requires columns ``total_counts``, ``n_genes_by_counts``, ``pct_counts_mt``.
    MAD bounds are computed on the count-floor-passing subset, then lower bounds
    are clamped to the absolute floors.

    Each ``*_override`` (when not ``None``) pins the *effective* filtering bound to
    that exact value; the MAD/floor-derived value is still computed and returned as
    the matching ``*_derived`` key. With no overrides the canonical keys equal the
    ``*_derived`` keys, so existing callers are unaffected.
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
        "total_counts_min": _eff(c_lo, total_counts_min_override),
        "total_counts_min_derived": float(c_lo),
        "total_counts_max": _eff(c_hi, total_counts_max_override),
        "total_counts_max_derived": float(c_hi),
        "total_counts_mad_lo_raw": c_lo_mad_raw,
        "n_genes_min": _eff(g_lo, n_genes_min_override),
        "n_genes_min_derived": float(g_lo),
        "n_genes_max": _eff(g_hi, n_genes_max_override),
        "n_genes_max_derived": float(g_hi),
        "n_genes_mad_lo_raw": g_lo_mad_raw,
        "pct_counts_mt_max": _eff(pct_mt_upper, pct_counts_mt_max_override),
        "pct_counts_mt_max_derived": float(pct_mt_upper),
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
    n_frag: np.ndarray,
    *,
    k_mad: float,
    n_frag_floor: float,
    n_fragments_min_override: float | None = None,
    n_fragments_max_override: float | None = None,
) -> tuple[float, float, float, tuple[float, float]]:
    """MAD bounds on log(n_fragments) after the absolute floor (mirrors S2).

    Returns ``(applied_lower, applied_upper, mad_lower_raw, (lower_derived,
    upper_derived))`` where ``mad_lower_raw`` is the log-MAD lower bound before the
    absolute ``n_frag_floor`` clamp, and the trailing pair is the MAD/floor-derived
    (pre-override) bounds. A user ``*_override`` pins the applied bound; the derived
    pair is unaffected so the figures can still draw the MAD line in grey. With no
    overrides ``applied_* == *_derived``.
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
    f_lo_derived, f_hi_derived = float(f_lo), float(f_hi)
    applied_lo = _eff(f_lo_derived, n_fragments_min_override)
    applied_hi = _eff(f_hi_derived, n_fragments_max_override)
    return applied_lo, applied_hi, f_lo_mad_raw, (f_lo_derived, f_hi_derived)


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
