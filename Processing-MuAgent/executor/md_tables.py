"""Shared markdown-table helpers.

Relocated out of ``qc_summary.py`` so both the plan-review / QC-exploration tables
and the post-QC summary render through one implementation. Formatting is identical
to the previous private helpers (``_fmt``, ``_md_table_cell``, ``_md_table``,
``_fmt_threshold_range``).
"""
from __future__ import annotations

from typing import Any

import numpy as np


def fmt(value: Any) -> str:
    """Format a scalar for a user-facing threshold table."""
    if value is None:
        return "n/a"
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (int, np.integer)):
        return f"{int(value)}"
    if isinstance(value, (float, np.floating)):
        v = float(value)
        if np.isnan(v):
            return "nan"
        if v == int(v) and abs(v) < 1e6:
            return f"{int(v)}"
        return f"{v:.2f}"
    if isinstance(value, (list, tuple)):
        return ", ".join(fmt(x) for x in value)
    return str(value)


def md_table_cell(value: Any) -> str:
    """Format a table cell; preserve inline markdown (**, trailing *)."""
    if isinstance(value, str) and ("**" in value or value.endswith("*")):
        return value
    return fmt(value)


def md_table(header: list[str], rows: list[list[Any]]) -> str:
    align = "|" + "|".join("---" for _ in header) + "|"
    h = "| " + " | ".join(md_table_cell(c) for c in header) + " |"
    body = "\n".join(
        "| " + " | ".join(md_table_cell(x) for x in r) + " |" for r in rows
    )
    return f"{h}\n{align}\n{body}"


def fmt_threshold_range(lo: Any, hi: Any) -> str:
    return f"{fmt(lo)} – {fmt(hi)}"
