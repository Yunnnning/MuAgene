"""MAD-based threshold helpers used in RNA / ATAC QC proposals."""
from __future__ import annotations

import numpy as np


def median_mad(values: np.ndarray) -> tuple[float, float]:
    med = float(np.median(values))
    mad = float(np.median(np.abs(values - med)))
    return med, mad


def upper_bound(values: np.ndarray, k: float = 3.0,
                floor: float | None = None, ceiling: float | None = None) -> float:
    med, mad = median_mad(np.asarray(values, dtype=float))
    b = med + k * mad
    if floor is not None:
        b = max(b, floor)
    if ceiling is not None:
        b = min(b, ceiling)
    return float(b)


def log_mad_bounds(values: np.ndarray, k: float = 5.0) -> tuple[float, float]:
    """Symmetric MAD bounds on log1p-transformed values; return linear lower/upper."""
    v = np.asarray(values, dtype=float)
    v = v[v > 0]
    lv = np.log1p(v)
    med, mad = median_mad(lv)
    lower = np.expm1(max(med - k * mad, 0.0))
    upper = np.expm1(med + k * mad)
    return float(lower), float(upper)
