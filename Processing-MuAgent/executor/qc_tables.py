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

# Upper-bound skip sentinels for 10x single-cell data.
# A threshold at or above these values is rendered as a skip label instead of
# the raw number — only the display changes, the filtering logic is unaffected.
# Public so qc_explore.py can import them for figure-line skip detection.
SKIP_ABOVE_COUNTS = 1_000_000  # total_counts / n_fragments: k_mad=999 yields values >> this
SKIP_ABOVE_GENES = 100_000     # n_genes: similar
SKIP_ABOVE_TSS = 500           # TSS enrichment (realistic max ~50–100)
SKIP_ABOVE_NUC = 50            # nucleosome signal ratio (realistic max ~3–5)
SKIP_PCT_AT = 100.0            # percent metrics: ≥ 100 means all cells pass


def _fmt_range_removal(lo, hi, skip_above: float) -> str:
    """Express a keep-if-in-range threshold as a removal condition.

    A cell is removed when its value falls outside [lo, hi].  With workaround
    skip values (hi >= skip_above) the upper bound is effectively infinite.

    Examples:
      lo=1000, hi=251464 → "< 1000 or > 251464"
      lo=0,    hi=251464 → "> 251464"
      lo=1000, hi=skip   → "< 1000"
      lo=0,    hi=skip   → "not applied"
    """
    lo_v = float(lo) if lo is not None else 0.0
    skipped = hi is not None and float(hi) >= skip_above
    if skipped:
        return "not applied" if lo_v <= 0.0 else f"< {_fmt(lo)}"
    if lo_v <= 0.0:
        return f"> {_fmt(hi)}"
    return f"< {_fmt(lo)} or > {_fmt(hi)}"


def _fmt_upper_removal(val, skip_above: float, prefix: str) -> str:
    """Express a keep-if-below threshold as a removal condition.

    Returns "not applied" when val >= skip_above (workaround skip sentinel),
    otherwise "{prefix}{fmt(val)}" (e.g. "> 5" or "≥ 3").
    """
    if val is not None and float(val) >= skip_above:
        return "not applied"
    return f"{prefix}{_fmt(val)}"


def _fmt_lower_removal(val, skip_at_or_below: float = 0.0) -> str:
    """Express a keep-if-above threshold as a removal condition.

    Returns "not applied" when val <= skip_at_or_below (user disabled the filter),
    otherwise "< {fmt(val)}".
    """
    if val is None or float(val) <= skip_at_or_below:
        return "not applied"
    return f"< {_fmt(val)}"


def frip_removal_condition(
    frip_min,
    *,
    peak_source: str | None = None,
    runtime_note: bool = False,
) -> str:
    """FRiP removal-condition label for QC / plan tables."""
    base = _fmt_lower_removal(frip_min)
    if base == "not applied":
        return base
    if runtime_note:
        return f"{base} _(computed at runtime)_"
    if not peak_source:
        return f"{base} _(not applied — no peaks available)_"
    return base


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
    value_label: str = "removed if",
    include_note: bool = False,
) -> str:
    """Render the RNA QC removal table.

    ``thresholds`` uses canonical keys ``total_counts_min/max``,
    ``n_genes_min/max``, ``pct_counts_mt_max``, ``pct_counts_ribo_max``.
    ``removals`` is a ``marginal_removals`` dict.
    """
    th = thresholds
    rows = [
        _row("total_counts",
             _fmt_range_removal(th.get("total_counts_min"), th.get("total_counts_max"), SKIP_ABOVE_COUNTS),
             _rm(removals, "total_counts"), include_note),
        _row("n_genes",
             _fmt_range_removal(th.get("n_genes_min"), th.get("n_genes_max"), SKIP_ABOVE_GENES),
             _rm(removals, "n_genes"), include_note),
        _row("pct_counts_mt",
             _fmt_upper_removal(th.get("pct_counts_mt_max"), SKIP_PCT_AT, "> "),
             _rm(removals, "pct_counts_mt"), include_note),
        _row("pct_counts_ribo",
             _fmt_upper_removal(th.get("pct_counts_ribo_max"), SKIP_PCT_AT, "> "),
             _rm(removals, "pct_counts_ribo"), include_note),
        _row("multiple_metrics", "—", _rm(removals, "multiple_metrics"), include_note),
        _row("total_removed", "—", _rm(removals, "total_removed"), include_note),
    ]
    return _render(value_label, rows, include_note)


def atac_removal_table(
    thresholds: dict[str, Any],
    removals: dict[str, Any],
    *,
    value_label: str = "removed if",
    include_note: bool = False,
    frip_threshold_display: str | None = None,
    frip_removed: Any = "",
    peak_source: str | None = None,
    frip_runtime_note: bool = False,
) -> str:
    """Render the ATAC QC removal table.

    ``thresholds`` uses canonical keys ``n_fragments_min/max``,
    ``tss_enrichment_min/max``, ``nucleosome_signal_max``, ``frip_min``.
    ``frip_threshold_display`` overrides the FRiP removal-condition cell (e.g.
    "< 0.20 _(computed at runtime)_"); ``frip_removed`` is the FRiP removal
    count (or ``"—"`` when FRiP was not applied). When ``frip_threshold_display``
    is omitted, ``frip_removal_condition`` derives the label from ``frip_min``,
    honoring user skip (``frip_min=0`` → "not applied") and optional peak /
    runtime annotations.
    """
    th = thresholds
    if frip_threshold_display is None:
        frip_threshold_display = frip_removal_condition(
            th.get("frip_min"),
            peak_source=peak_source,
            runtime_note=frip_runtime_note,
        )
    rows = [
        _row("n_fragments",
             _fmt_range_removal(th.get("n_fragments_min"), th.get("n_fragments_max"), SKIP_ABOVE_COUNTS),
             _rm(removals, "n_fragments"), include_note),
        _row("tss_enrichment",
             _fmt_range_removal(th.get("tss_enrichment_min"), th.get("tss_enrichment_max"),
                                SKIP_ABOVE_TSS),
             _rm(removals, "tss_enrichment"), include_note),
        _row("nucleosome_signal",
             _fmt_upper_removal(th.get("nucleosome_signal_max"), SKIP_ABOVE_NUC, "≥ "),
             _rm(removals, "nucleosome_signal"), include_note),
        _row("frip", frip_threshold_display, frip_removed, include_note),
        _row("multiple_metrics", "—", _rm(removals, "multiple_metrics"), include_note),
        _row("total_removed", "—", _rm(removals, "total_removed"), include_note),
    ]
    return _render(value_label, rows, include_note)
