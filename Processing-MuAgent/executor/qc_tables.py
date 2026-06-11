"""Shared QC removal-table builders.

One implementation of the RNA / ATAC "cells removed" tables, used by both the
pre-plan QC exploration appendix (``executor.qc_explore``) and the post-filter
QC review/summary (``executor.qc_summary``). Keeping a single builder guarantees
both reports present the same marginal removal counts with identical layout.

Counts are *marginal* (see ``methods.qc_filter_stats.marginal_removals``): each
per-metric row counts every cell failing that one threshold; ``multiple_metrics``
counts cells failing two or more; ``total_removed`` is the union.
"""
from __future__ import annotations

from typing import Any

from .md_tables import fmt as _fmt
from .md_tables import fmt_threshold_range as _fmt_range
from .md_tables import md_table as _md_table

# One-line, user-facing explanations of each metric / shortcode for the optional
# "note" column. Keys match the parameter-column values used below.
NOTES: dict[str, str] = {
    "total_counts": "Total UMI counts per cell (library size); low = empty/dying, high = potential doublets.",
    "n_genes": "Genes detected per cell; low = low-quality, high = potential doublets.",
    "pct_counts_mt": "Percent of counts from mitochondrial genes — high indicates stressed or dying cells.",
    "pct_counts_ribo": "Percent of counts from ribosomal protein genes (Rps/Rpl/Mrps/Mrpl).",
    "n_fragments": "ATAC fragments per cell (library depth).",
    "tss_enrichment": "Fragment enrichment at transcription start sites (signal-to-noise).",
    "nucleosome_signal": "Mono- to nucleosome-free fragment ratio (nucleosome positioning quality).",
    "frip": "Fraction of reads in peaks — computed at runtime when a peak set is available.",
    "multiple_metrics": "Cells failing two or more thresholds (counted once).",
    "total_removed": "Cells removed by the combined filter (union of all thresholds).",
}


def _rm(removals: dict[str, Any], key: str) -> Any:
    v = removals.get(key)
    return v if v is not None else ""


def _row(param: str, value: Any, removed: Any, include_note: bool) -> list[Any]:
    row = [param, value, removed]
    if include_note:
        row.append(NOTES.get(param, ""))
    return row


def _render(value_label: str, rows: list[list[Any]], include_note: bool) -> str:
    header = ["parameter", value_label, "cells removed"]
    if include_note:
        header.append("note")
    return _md_table(header, rows)


def rna_removal_table(
    thresholds: dict[str, Any],
    removals: dict[str, Any],
    *,
    value_label: str = "value",
    include_note: bool = False,
) -> str:
    """Render the RNA QC removal table.

    ``thresholds`` uses canonical keys ``total_counts_min/max``,
    ``n_genes_min/max``, ``pct_counts_mt_max``, ``pct_counts_ribo_max``.
    ``removals`` is a ``marginal_removals`` dict.
    """
    th = thresholds
    rows = [
        _row("total_counts", _fmt_range(th.get("total_counts_min"), th.get("total_counts_max")),
             _rm(removals, "total_counts"), include_note),
        _row("n_genes", _fmt_range(th.get("n_genes_min"), th.get("n_genes_max")),
             _rm(removals, "n_genes"), include_note),
        _row("pct_counts_mt", f"≤ {_fmt(th.get('pct_counts_mt_max'))}",
             _rm(removals, "pct_counts_mt"), include_note),
        _row("pct_counts_ribo", f"≤ {_fmt(th.get('pct_counts_ribo_max'))}",
             _rm(removals, "pct_counts_ribo"), include_note),
        _row("multiple_metrics", "—", _rm(removals, "multiple_metrics"), include_note),
        _row("total_removed", "—", _rm(removals, "total_removed"), include_note),
    ]
    return _render(value_label, rows, include_note)


def atac_removal_table(
    thresholds: dict[str, Any],
    removals: dict[str, Any],
    *,
    value_label: str = "value",
    include_note: bool = False,
    frip_threshold_display: str | None = None,
    frip_removed: Any = "",
) -> str:
    """Render the ATAC QC removal table.

    ``thresholds`` uses canonical keys ``n_fragments_min/max``,
    ``tss_enrichment_min/max``, ``nucleosome_signal_max``, ``frip_min``.
    ``frip_threshold_display`` overrides the FRiP value cell (e.g. a "computed at
    runtime" or "not applied — no peaks" annotation); ``frip_removed`` is the FRiP
    removal count (or ``"—"`` when FRiP was not applied).
    """
    th = thresholds
    if frip_threshold_display is None:
        frip_threshold_display = f"≥ {_fmt(th.get('frip_min'))}"
    rows = [
        _row("n_fragments", _fmt_range(th.get("n_fragments_min"), th.get("n_fragments_max")),
             _rm(removals, "n_fragments"), include_note),
        _row("tss_enrichment", _fmt_range(th.get("tss_enrichment_min"), th.get("tss_enrichment_max")),
             _rm(removals, "tss_enrichment"), include_note),
        _row("nucleosome_signal", f"< {_fmt(th.get('nucleosome_signal_max'))}",
             _rm(removals, "nucleosome_signal"), include_note),
        _row("frip", frip_threshold_display, frip_removed, include_note),
        _row("multiple_metrics", "—", _rm(removals, "multiple_metrics"), include_note),
        _row("total_removed", "—", _rm(removals, "total_removed"), include_note),
    ]
    return _render(value_label, rows, include_note)
