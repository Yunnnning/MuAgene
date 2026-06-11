"""Order-independent marginal QC removal counts.

Each metric row counts *every* cell that fails that one threshold, evaluated
independently against the full unfiltered set. The same cell can therefore be
counted under several metrics (intended). ``multiple_metrics`` counts cells
failing two or more metrics; ``total_removed`` is the union (cells failing at
least one). Counts are order-independent.

This is the single counting standard shared by the pre-plan QC exploration
(``executor.qc_explore``) and the real QC stages (``s1_rna_qc`` / ``s2_atac_qc``)
so the planning preview and the post-filter report use the same semantics.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def marginal_removals(pass_masks: dict[str, np.ndarray]) -> dict[str, int]:
    """Return non-exclusive per-metric removal counts.

    Each metric counts *every* cell that fails that one threshold, evaluated
    independently against the full array. The same cell can be counted under
    several metrics (intended). Two summary keys are added:

    - ``multiple_metrics``: cells failing two or more metrics.
    - ``total_removed``: cells failing at least one metric (the AND-filter union).

    Counts are order-independent: each mask is applied to the unfiltered set.
    """
    if not pass_masks:
        return {"multiple_metrics": 0, "total_removed": 0}

    names = list(pass_masks.keys())
    n = int(next(iter(pass_masks.values())).size)
    fail_count = np.zeros(n, dtype=int)
    out: dict[str, int] = {}
    for name in names:
        fail_m = ~np.asarray(pass_masks[name], dtype=bool)
        out[name] = int(fail_m.sum())
        fail_count += fail_m.astype(int)

    out["multiple_metrics"] = int((fail_count >= 2).sum())
    out["total_removed"] = int((fail_count >= 1).sum())
    return out


def append_frip(
    counts: dict[str, Any],
    *,
    frip_fail: int | None,
    n_pre: int,
    n_post: int,
) -> dict[str, Any]:
    """Add FRiP-only removals and refresh total_removed for the full S2 filter.

    FRiP is evaluated only on cells that pass the three core metrics, so FRiP
    failures are disjoint from the core-metric failures and do not change
    ``multiple_metrics``. ``total_removed`` becomes the overall union
    ``n_pre - n_post`` (core union plus FRiP failures).
    """
    result = dict(counts)
    if frip_fail is not None:
        result["frip_min"] = int(frip_fail)
    result["total_removed"] = int(n_pre - n_post)
    return result
