"""Order-independent exclusive QC removal counts.

Each metric row counts cells that fail only that metric while passing all
others in the group. Cells failing multiple metrics are counted under
``multiple_metrics``. ``total_removed`` is the combined AND-filter removal count.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def exclusive_removals(pass_masks: dict[str, np.ndarray]) -> dict[str, int]:
    """Return per-metric exclusive counts plus multiple_metrics and total_removed."""
    if not pass_masks:
        return {"multiple_metrics": 0, "total_removed": 0}

    names = list(pass_masks.keys())
    n = int(next(iter(pass_masks.values())).size)
    combined = np.ones(n, dtype=bool)
    for pass_m in pass_masks.values():
        combined &= np.asarray(pass_m, dtype=bool)

    out: dict[str, int] = {}
    for name in names:
        pass_m = np.asarray(pass_masks[name], dtype=bool)
        others = np.ones(n, dtype=bool)
        for other_name, other_pass in pass_masks.items():
            if other_name != name:
                others &= np.asarray(other_pass, dtype=bool)
        out[name] = int((~pass_m & others).sum())

    n_fail = int((~combined).sum())
    out["multiple_metrics"] = n_fail - sum(out[n] for n in names)
    out["total_removed"] = n_fail
    return out


def append_frip_exclusive(
    counts: dict[str, Any],
    *,
    frip_fail: int | None,
    n_pre: int,
    n_post: int,
) -> dict[str, Any]:
    """Add FRiP-only removals and refresh total_removed for the full S2 filter."""
    result = dict(counts)
    if frip_fail is not None:
        result["frip_min"] = int(frip_fail)
    result["total_removed"] = int(n_pre - n_post)
    return result
