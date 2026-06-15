"""Figure helpers: consistent size, font, DPI; PNG + PDF (vector) companion.

Design goals:
  - PNG at >= 300 dpi for raster viewing.
  - PDF (vector) companion written at the same path with .pdf extension.
  - Larger readable fonts; fixed figure size per figure type for consistency.
  - Deterministic, stage-prefixed filenames.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


FIGURE_DPI = 300
FONT_SIZE = 14
QC_FILL_COLOR = "#f97316"
QC_FILL_ALPHA = 0.45
QC_EDGE_COLOR = "#c2410c"
QC_EDGE_LINEWIDTH = 0.4
RNA_HIST_FILL_COLOR = "royalblue"
RNA_HIST_FILL_ALPHA = 0.50
RNA_HIST_EDGE_COLOR = "#1e40af"
QC_ANNOTATION_COLOR = "firebrick"
ANNOTATION_LINEWIDTH = 1.2
ANNOTATION_FONTSIZE = FONT_SIZE - 1
TITLE_SIZE = 16
RNA_VIOLIN_FILL_COLOR = "dodgerblue"
RNA_VIOLIN_FILL_ALPHA = 0.50
RNA_VIOLIN_QUANTILE_COLOR = "mediumblue"
RNA_VIOLIN_QUANTILE_LINEWIDTH = 1.6
FRIP_DISTRIBUTION_TITLE = "Fraction of Reads in Peaks (FRiP) distribution"
FRIP_XMAX = 0.7
QC_PLOT_PAIR_SIZE = (7.0, 4.5)
TSS_PROFILE_TITLE = "Mean TSS enrichment scores"
TSS_PROFILE_CAPTION = (
    "left = cells passing the TSS enrichment threshold, "
    "right = cells failing the TSS enrichment threshold."
)
TSS_PROFILE_WINDOW_BP = 1500
TSS_PASS_COLOR = "coral"
TSS_FAIL_COLOR = "lightseagreen"
UMAP_SIZE = (6.5, 5.5)
QC_VIOLIN_SIZE = (12, 4.5)
QC_REF_LINE_COLOR = "#6b7280"  # grey for fixed reference markers (e.g. mt ceilings)
THRESHOLD_LABEL_Y_TOP = 0.98
THRESHOLD_LABEL_Y_STEP = 0.13
THRESHOLD_LABEL_Y_MIN = 0.32
THRESHOLD_LABEL_X_SEP_FRAC = 0.05  # fraction of x-range treated as "same slot"
THRESHOLD_LABEL_X_OFFSET_LINEAR_FRAC = 0.018  # label shift off vline (linear axes)
THRESHOLD_LABEL_X_OFFSET_PCT = 0.35  # fixed shift for % metrics (robust to outlier x-range)
THRESHOLD_LABEL_X_OFFSET_LOG_FRAC = 0.055  # label shift off vline (log axes)
QC_EXPLORE_RNA_TITLE = "RNA Data Distribution and QC Thresholds"
QC_EXPLORE_ATAC_TITLE = "ATAC Data Distribution and QC Thresholds"
QC_HIST_SUPTITLE_SIZE = TITLE_SIZE + 1
QC_PANEL_SUBTITLE_SIZE = FONT_SIZE
QC_HIST_METRIC_TITLE_Y = 1.12  # transAxes; va=bottom
QC_HIST_REMOVED_SUBTITLE_Y = 1.04  # transAxes; va=bottom — gap above plot top (1.0)
QC_HIST_Y_TOP_PAD = 1.10  # expand ylim after hist so bars clear the subtitle band
QC_HIST_PANEL_W = 5.5
QC_HIST_PANEL_H = 4.6
# The per-panel metric name + "(cells removed: N)" are floating ax.text at
# transAxes y>1.0 (see _set_qc_panel_titles), which tight_layout does NOT account
# for. Reserve a generous top band so they clear the suptitle, and add row spacing
# so the bottom row's floating titles clear the top row's x-axis labels.
QC_HIST_SUPTITLE_Y = 0.975
QC_HIST_LAYOUT_RECT = (0, 0, 1, 0.863)
QC_HIST_SUBPLOTS_TOP = 0.863
QC_HIST_HSPACE = 0.50

# User-facing panel / axis labels for QC explore histograms (keys = internal metric ids).
QC_METRIC_DISPLAY_NAMES: dict[str, str] = {
    "total_counts": "total counts",
    "n_genes": "number of genes",
    "pct_counts_mt": "mitochondrial percentage",
    "pct_counts_ribo": "ribosomal percentage",
    "n_fragments": "number of fragments",
    "tss_enrichment": "TSS enrichment",
    "nucleosome_signal": "nucleosome signal",
}


def _qc_metric_display_name(name: str) -> str:
    return QC_METRIC_DISPLAY_NAMES.get(name, name.replace("_", " "))


def _qc_hist_grid_shape(n: int) -> tuple[int, int]:
    """Return (nrows, ncols) for QC threshold histogram panels."""
    if n <= 1:
        return 1, 1
    if n == 2:
        return 1, 2
    ncols = 2
    return (n + ncols - 1) // ncols, ncols


def _set_qc_panel_titles(ax, name: str, n_removed: int) -> None:
    """Metric name + smaller removed-count subtitle above the histogram."""
    ax.set_title("")
    ax.text(
        0.5, QC_HIST_METRIC_TITLE_Y, name,
        transform=ax.transAxes, ha="center", va="bottom",
        fontsize=TITLE_SIZE, clip_on=False,
    )
    ax.text(
        0.5, QC_HIST_REMOVED_SUBTITLE_Y, f"(cells removed: {n_removed:,})",
        transform=ax.transAxes, ha="center", va="bottom",
        fontsize=QC_PANEL_SUBTITLE_SIZE, color="#333333", clip_on=False,
    )


def _is_pct_metric(name: str) -> bool:
    return "pct" in name.lower()


def _thresholds_coincide(a: float, b: float | None, *, pct: bool) -> bool:
    """True when two cutoffs are the same marker (skip duplicate ref line)."""
    if b is None:
        return False
    tol = 0.2 if pct else max(1e-9, abs(b) * 0.005)
    return abs(a - b) <= tol


def _format_cutoff_value(x: float, *, pct: bool, log_axis: bool) -> str:
    if pct:
        rounded = round(x)
        if abs(x - rounded) < 0.05:
            return f"{int(rounded)}%"
        return f"{x:.1f}%"
    if log_axis:
        return f"{x:,.0f}" if x >= 100 else f"{x:.2g}"
    return f"{x:.2g}"


def _cutoff_label(
    x: float, *, pct: bool, log_axis: bool, mad: bool = False,
) -> str:
    text = _format_cutoff_value(x, pct=pct, log_axis=log_axis)
    return f"{text} (MAD)" if mad else text


def _stagger_threshold_label_ys(
    xs: list[float],
    x_range: float,
    *,
    active: list[bool] | None = None,
) -> list[float]:
    """Assign axes-fraction y anchors so rotated labels on nearby cutoffs don't overlap.

    Chosen (active) thresholds always anchor at the top; only reference markers
    are staggered when their x positions cluster.
    """
    n = len(xs)
    if n == 0:
        return []
    if active is None:
        active = [False] * n

    ys = [THRESHOLD_LABEL_Y_TOP] * n
    inactive = [i for i, is_active in enumerate(active) if not is_active]
    if not inactive:
        return ys

    min_sep = max(x_range * THRESHOLD_LABEL_X_SEP_FRAC, 1e-12)
    order = sorted(inactive, key=lambda i: xs[i])
    cluster_slot = 0
    for rank, idx in enumerate(order):
        if rank > 0 and abs(xs[idx] - xs[order[rank - 1]]) < min_sep:
            cluster_slot += 1
        else:
            cluster_slot = 0
        ys[idx] = max(
            THRESHOLD_LABEL_Y_MIN,
            THRESHOLD_LABEL_Y_TOP - cluster_slot * THRESHOLD_LABEL_Y_STEP,
        )
    return ys


def _label_x_offset(x: float, x_range: float, *, log_axis: bool, pct: bool = False) -> float:
    """Shift label anchor off the vline so the line does not bisect the text."""
    if log_axis and x > 0:
        return x * (1.0 + THRESHOLD_LABEL_X_OFFSET_LOG_FRAC)
    if pct:
        return x + THRESHOLD_LABEL_X_OFFSET_PCT
    return x + max(x_range * THRESHOLD_LABEL_X_OFFSET_LINEAR_FRAC, 1e-12)


def _draw_threshold_markers(
    ax,
    markers: list[tuple[float, str, bool]],
    *,
    x_range: float,
    log_axis: bool = False,
    pct: bool = False,
) -> None:
    """Draw cutoff vlines and staggered vertical labels.

    Each marker is ``(x, label, is_active)``; active lines/labels use
    ``QC_ANNOTATION_COLOR``, reference markers use ``QC_REF_LINE_COLOR``.
    """
    if not markers:
        return
    xs = [m[0] for m in markers]
    actives = [m[2] for m in markers]
    ys = _stagger_threshold_label_ys(xs, x_range, active=actives)
    for (x, label, is_active), y_frac in zip(markers, ys):
        color = QC_ANNOTATION_COLOR if is_active else QC_REF_LINE_COLOR
        linestyle = "--" if is_active else ":"
        ax.axvline(
            x, color=color, linestyle=linestyle,
            linewidth=ANNOTATION_LINEWIDTH, zorder=5 if is_active else 4,
        )
        label_x = _label_x_offset(x, x_range, log_axis=log_axis, pct=pct)
        ax.text(
            label_x, y_frac, f" {label}",
            transform=ax.get_xaxis_transform(),
            rotation=90, va="top", ha="left",
            fontsize=ANNOTATION_FONTSIZE - 1, color=color,
            zorder=6, clip_on=False,
        )


def _apply_style():
    import matplotlib as mpl
    mpl.use("Agg")
    mpl.rcParams.update({
        "font.size": FONT_SIZE,
        "axes.titlesize": TITLE_SIZE,
        "axes.labelsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE - 1,
        "ytick.labelsize": FONT_SIZE - 1,
        "legend.fontsize": FONT_SIZE - 1,
        "savefig.dpi": FIGURE_DPI,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,  # TrueType so text is editable in the PDF
        "svg.fonttype": "none",
    })


def save_figure(
    fig,
    out_dir: Path | str,
    stem: str,
    *,
    also_pdf: bool = True,
    dpi: int | None = None,
) -> list[Path]:
    """Save `fig` as <stem>.png and optionally <stem>.pdf. Returns list of paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    png_path = out_dir / f"{stem}.png"
    if not also_pdf:
        pdf_stale = out_dir / f"{stem}.pdf"
        if pdf_stale.exists():
            pdf_stale.unlink()
    fig.savefig(png_path, dpi=dpi or FIGURE_DPI, bbox_inches="tight")
    paths.append(png_path)
    if also_pdf:
        pdf_path = out_dir / f"{stem}.pdf"
        fig.savefig(pdf_path, bbox_inches="tight")
        paths.append(pdf_path)
    import matplotlib.pyplot as plt
    plt.close(fig)
    from . import io as _io
    for path in paths:
        _io.sync_path(path)
    return paths


def plot_umap(coords: np.ndarray, labels, *, title: str, out_dir: Path | str,
              stem: str, label_name: str = "cluster") -> list[Path]:
    """Render a 2D UMAP scatter coloured by discrete labels."""
    _apply_style()
    import matplotlib.pyplot as plt
    coords = np.asarray(coords)
    labels_arr = np.asarray(labels)
    fig, ax = plt.subplots(figsize=UMAP_SIZE)
    uniq = sorted(set(map(str, labels_arr)), key=lambda s: (len(s), s))
    cmap = plt.get_cmap("tab20" if len(uniq) > 10 else "tab10")
    for i, v in enumerate(uniq):
        mask = labels_arr.astype(str) == v
        ax.scatter(coords[mask, 0], coords[mask, 1], s=8,
                   color=cmap(i % cmap.N), label=str(v), alpha=0.85, linewidths=0)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(title=label_name, fontsize=FONT_SIZE - 2, markerscale=1.8,
              bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return save_figure(fig, out_dir, stem)


def plot_qc_violin(values_dict: dict[str, np.ndarray], *, out_dir: Path | str,
                    stem: str, title: str) -> list[Path]:
    """Render side-by-side violins for 1..N QC metrics."""
    _apply_style()
    import matplotlib.pyplot as plt
    n = len(values_dict)
    # Scale figure width with the metric count so violins stay readable.
    fig_w = max(QC_VIOLIN_SIZE[0], 3.0 * n)
    fig, axes = plt.subplots(1, n, figsize=(fig_w, QC_VIOLIN_SIZE[1]))
    if n == 1:
        axes = [axes]
    for ax, (name, vals) in zip(axes, values_dict.items()):
        v = np.asarray(vals, dtype=float)
        v = v[np.isfinite(v)]
        label = _qc_metric_display_name(name)
        if v.size == 0:
            ax.set_title(label + " (no data)")
            continue
        parts = ax.violinplot(v, showmeans=False, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(RNA_VIOLIN_FILL_COLOR)
            pc.set_edgecolor(RNA_VIOLIN_QUANTILE_COLOR)
            pc.set_alpha(RNA_VIOLIN_FILL_ALPHA)
        for key in ("cmedians", "cbars", "cmins", "cmaxes"):
            if key in parts:
                parts[key].set_color(RNA_VIOLIN_QUANTILE_COLOR)
                parts[key].set_linewidth(RNA_VIOLIN_QUANTILE_LINEWIDTH)
        ax.set_title(label)
        ax.set_ylabel(label)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(title)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return save_figure(fig, out_dir, stem)


def plot_qc_threshold_histograms(
    metrics: "dict[str, dict[str, Any]]",
    *,
    out_dir: Path | str,
    stem: str,
    title: str,
    fill_color: str = QC_FILL_COLOR,
    edge_color: str = QC_EDGE_COLOR,
    fill_alpha: float = QC_FILL_ALPHA,
    extra_panel: "dict[str, Any] | None" = None,
    extra_panel_slot: "int | None" = None,
) -> list[Path]:
    """Per-metric histograms with the default QC cutoffs drawn on them.

    ``metrics`` maps a metric name to a spec dict:
      - ``values``: 1D array of per-cell values.
      - ``lo`` / ``hi``: MAD-derived cutoffs (either may be ``None``); cells
        outside ``[lo, hi]`` are shaded as removed and the marginal count is
        annotated.
      - ``log``: when truthy, use a log-spaced x-axis (counts / fragments).
      - ``refs``: optional list of ``(x, label)`` fixed reference lines drawn in a
        distinct style (e.g. 5% / 10% mito ceilings) — markers only, they do not
        affect the removed count.
      - ``mad_lo_raw`` / ``mad_hi_raw``: optional raw MAD bounds; active ``lo`` /
        ``hi`` labels get a `` (MAD)`` suffix when they coincide with the raw value.
      - ``mad_hi``: when true and ``mad_hi_raw`` is omitted, the upper cutoff is
        always labeled as MAD-derived (log-MAD upper bounds without clamping).

    ``extra_panel`` (optional) fills one grid slot with a non-histogram panel:
    ``{"distr": <frag_size_distr vector>, "title": str}``. ``extra_panel_slot``
    is the 0-based slot index (row-major in the 2×2 ATAC grid); defaults to the
    last slot when omitted.
    """
    _apply_style()
    import matplotlib.pyplot as plt

    n = len(metrics)
    if n == 0:
        fig, ax = plt.subplots(figsize=QC_PLOT_PAIR_SIZE)
        ax.set_title(title + " (no data)")
        return save_figure(fig, out_dir, stem)

    n_used = n + (1 if extra_panel else 0)
    nrows, ncols = _qc_hist_grid_shape(n_used)
    fig_w = QC_HIST_PANEL_W * ncols
    fig_h = QC_HIST_PANEL_H * nrows
    fig, axes_grid = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h))
    axes = list(np.atleast_1d(axes_grid).ravel())
    for ax in axes[n_used:]:
        ax.set_visible(False)

    extra_slot = extra_panel_slot if extra_panel_slot is not None else n
    if extra_panel is not None:
        if extra_slot < 0 or extra_slot >= n_used:
            raise ValueError(
                f"extra_panel_slot={extra_slot} out of range for {n_used} panels"
            )
        hist_slots = [i for i in range(n_used) if i != extra_slot]
        _draw_fragment_size_distribution(
            axes[extra_slot], extra_panel.get("distr"),
            title=extra_panel.get("title"),
        )
    else:
        hist_slots = list(range(n))

    for slot, (name, spec) in zip(hist_slots, metrics.items()):
        ax = axes[slot]
        vals = np.asarray(spec.get("values"), dtype=float)
        vals = vals[np.isfinite(vals)]
        lo = spec.get("lo")
        hi = spec.get("hi")
        refs = spec.get("refs") or []
        use_log = bool(spec.get("log"))
        if vals.size == 0:
            ax.set_title(_qc_metric_display_name(name) + " (no data)")
            continue

        if use_log:
            pos = vals[vals > 0]
            if pos.size:
                bins = np.logspace(np.log10(pos.min()), np.log10(pos.max()), 50)
                ax.set_xscale("log")
            else:
                bins = 50
        else:
            bins = 50
        ax.hist(vals, bins=bins, color=fill_color, alpha=fill_alpha,
                edgecolor=edge_color, linewidth=QC_EDGE_LINEWIDTH)
        ymin, ymax = ax.get_ylim()
        if ymax > 0:
            ax.set_ylim(ymin, ymax * QC_HIST_Y_TOP_PAD)

        pct_metric = _is_pct_metric(name)
        mad_lo_raw = spec.get("mad_lo_raw")
        mad_hi_raw = spec.get("mad_hi_raw")
        removed = np.zeros(vals.shape, dtype=bool)
        markers: list[tuple[float, str, bool]] = []
        if lo is not None:
            removed |= vals < lo
            lo_mad = (
                mad_lo_raw is not None
                and _thresholds_coincide(float(mad_lo_raw), float(lo), pct=pct_metric)
            )
            markers.append((
                float(lo),
                _cutoff_label(float(lo), pct=pct_metric, log_axis=use_log, mad=lo_mad),
                True,
            ))
        if hi is not None:
            removed |= vals > hi
            if mad_hi_raw is not None:
                hi_mad = _thresholds_coincide(float(mad_hi_raw), float(hi), pct=pct_metric)
            else:
                hi_mad = bool(spec.get("mad_hi", False))
            markers.append((
                float(hi),
                _cutoff_label(float(hi), pct=pct_metric, log_axis=use_log, mad=hi_mad),
                True,
            ))
        for rx, rlabel in refs:
            rx_f = float(rx)
            if any(_thresholds_coincide(rx_f, b, pct=pct_metric) for b in (lo, hi)):
                continue
            markers.append((rx_f, str(rlabel), False))

        if use_log and vals.size:
            pos = vals[vals > 0]
            x_range = (
                float(np.log10(pos.max()) - np.log10(pos.min()))
                if pos.size else 1.0
            )
        else:
            x_range = float(vals.max() - vals.min()) if vals.size else 1.0
        _draw_threshold_markers(
            ax, markers, x_range=max(x_range, 1e-12), log_axis=use_log,
            pct=pct_metric,
        )

        _set_qc_panel_titles(ax, _qc_metric_display_name(name), int(removed.sum()))
        ax.set_xlabel(_qc_metric_display_name(name))
        ax.set_ylabel("number of cells")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout(rect=QC_HIST_LAYOUT_RECT)
    fig.subplots_adjust(top=QC_HIST_SUBPLOTS_TOP, hspace=QC_HIST_HSPACE)
    fig.suptitle(title, fontsize=QC_HIST_SUPTITLE_SIZE, y=QC_HIST_SUPTITLE_Y)
    return save_figure(fig, out_dir, stem)


def plot_fragment_size_distribution(
    distr: np.ndarray,
    *,
    out_dir: Path | str,
    stem: str,
    title: str,
    distr_after: np.ndarray | None = None,
) -> list[Path]:
    """Plot a 1D fragment-size histogram (counts per fragment length).

    `distr` is a 1D vector where `distr[i]` is the number of fragments of length
    `i` (the SnapATAC2 `frag_size_distr` format; index 0 holds the over-max bucket).
    A well-prepared ATAC library shows distinct peaks at ~150, ~300, ~450 bp.

    When `distr_after` is supplied, only the post-QC distribution is plotted
    (normalised to fraction of fragments).
    """
    _apply_style()
    import matplotlib.pyplot as plt

    # Prefer the post-QC distribution when supplied (the standalone S2 figure
    # plots a single curve); otherwise the input distribution.
    chosen = distr
    if distr_after is not None and np.asarray(distr_after).size > 0:
        chosen = distr_after

    fig, ax = plt.subplots(figsize=QC_PLOT_PAIR_SIZE)
    _draw_fragment_size_distribution(ax, chosen, title=title)
    return save_figure(fig, out_dir, stem)


def _draw_fragment_size_distribution(
    ax, distr: np.ndarray, *, title: str | None = None, x_right_max: int = 1000,
) -> None:
    """Render a normalised fragment-size distribution onto ``ax``.

    Stateless rendering shared by the standalone S2 figure and the QC-explore
    grid's 4th panel — the two callers compute their own distributions (pre- vs
    post-QC cell sets) and must not share data.

    ``distr`` is the SnapATAC2 ``frag_size_distr`` vector (index 0 is the
    over-max bucket and is dropped).
    """
    d = np.asarray(distr, dtype=float).ravel()
    if d.size == 0:
        ax.set_title(((title + " ") if title else "") + "(no data)")
        return
    body = d[1:]  # drop over-max bucket at index 0
    total = body.sum()
    normed = body / total if total > 0 else body
    x_full = np.arange(1, body.size + 1)
    x_right = min(x_right_max, body.size)
    x = x_full[:x_right]
    y = normed[:x_right]

    ax.fill_between(x, 0, y, alpha=QC_FILL_ALPHA, color=QC_FILL_COLOR, linewidth=0)
    ax.plot(x, y, color=QC_EDGE_COLOR, linewidth=QC_EDGE_LINEWIDTH)
    ax.set_xlim(left=0, right=x_right)
    ax.set_ylabel("fraction of fragments")
    ax.set_xlabel("fragment length (bp)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if title:
        ax.set_title(title)


def plot_frip_histogram(
    frip: np.ndarray,
    *,
    out_dir: Path | str,
    stem: str,
    frip_min: float = 0.25,
) -> list[Path]:
    """Histogram of per-cell FRiP values with a vertical threshold line.

    `frip` contains FRiP values for all cells that passed the 3-metric QC
    (before the FRiP filter itself). The vertical line shows the chosen threshold.
    """
    _apply_style()
    import matplotlib.pyplot as plt

    frip = np.asarray(frip, dtype=float)
    frip = frip[np.isfinite(frip)]

    fig, ax = plt.subplots(figsize=QC_PLOT_PAIR_SIZE)
    if frip.size == 0:
        ax.set_title(f"{FRIP_DISTRIBUTION_TITLE} (no data)")
        return save_figure(fig, out_dir, stem)

    bins = np.linspace(0, FRIP_XMAX, 36)
    ax.hist(
        frip, bins=bins, color=QC_FILL_COLOR, alpha=QC_FILL_ALPHA,
        edgecolor=QC_EDGE_COLOR, linewidth=QC_EDGE_LINEWIDTH,
    )
    x_range = float(frip.max() - frip.min()) if frip.size else 1.0
    _draw_threshold_markers(
        ax,
        [(frip_min, _cutoff_label(frip_min, pct=False, log_axis=False), True)],
        x_range=max(x_range, 1e-12),
    )
    ax.set_xlabel("FRiP (fraction of reads in peaks)")
    ax.set_ylabel("number of cells")
    ax.set_title(FRIP_DISTRIBUTION_TITLE)
    ax.set_xlim(0, FRIP_XMAX)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return save_figure(fig, out_dir, stem)


def _normalize_tss_profile(raw_profile: np.ndarray) -> np.ndarray:
    """Fold enrichment over flanking background (±2 kb, excluding central 200 bp)."""
    prof = np.asarray(raw_profile, dtype=float).ravel()
    if prof.size == 0:
        return prof
    center = prof.size // 2
    flank = np.concatenate([prof[: center - 100], prof[center + 101 :]])
    denom = float(flank.mean()) if flank.size and flank.mean() > 0 else 1.0
    return prof / denom


def plot_tss_enrichment_profile(
    profile_pass: np.ndarray,
    profile_fail: np.ndarray,
    *,
    out_dir: Path | str,
    stem: str,
    n_pass: int | None = None,
    n_fail: int | None = None,
) -> list[Path]:
    """Side-by-side mean TSS enrichment profiles for cells passing vs failing the TSS enrichment threshold."""
    _apply_style()
    import matplotlib.pyplot as plt

    pass_prof = _normalize_tss_profile(profile_pass)
    fail_prof = _normalize_tss_profile(profile_fail)
    n_bins = min(pass_prof.size, fail_prof.size)
    if n_bins == 0:
        fig, ax = plt.subplots(figsize=(8.0, 4.5))
        ax.set_title(f"{TSS_PROFILE_TITLE} (no data)")
        return save_figure(fig, out_dir, stem)

    pass_prof = pass_prof[:n_bins]
    fail_prof = fail_prof[:n_bins]
    center = n_bins // 2
    win = min(TSS_PROFILE_WINDOW_BP, center)
    sl = slice(center - win, center + win + 1)
    pass_prof = pass_prof[sl]
    fail_prof = fail_prof[sl]
    x = np.arange(-win, win + 1)

    pass_label = (
        f"pass TSS threshold ({n_pass:,} cells)"
        if n_pass is not None else "pass TSS threshold"
    )
    fail_label = (
        f"fail TSS threshold ({n_fail:,} cells)"
        if n_fail is not None else "fail TSS threshold"
    )
    tss_line_alpha = 0.88
    tss_line_width = 1.1
    y_top = float(max(pass_prof.max(initial=0), fail_prof.max(initial=0), 1.0))

    fig_w = QC_PLOT_PAIR_SIZE[0] * 2 + 0.5
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, QC_PLOT_PAIR_SIZE[1]), sharey=True)
    panels = (
        (axes[0], pass_prof, TSS_PASS_COLOR, pass_label),
        (axes[1], fail_prof, TSS_FAIL_COLOR, fail_label),
    )
    for ax, prof, color, title in panels:
        ax.plot(
            x, prof, color=color, linewidth=tss_line_width,
            alpha=tss_line_alpha, zorder=3,
        )
        ax.axvline(0, color="black", linewidth=1.0, linestyle="--", alpha=0.6, zorder=1)
        ax.set_xlim(-win, win)
        ax.set_ylim(0, y_top * 1.05)
        ax.set_xlabel("distance from TSS (bp)")
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("mean TSS enrichment score")
    fig.suptitle(TSS_PROFILE_TITLE, y=1.02)
    fig.tight_layout()
    return save_figure(fig, out_dir, stem)


MARKER_TSNE_MAIN_TITLE = "Marker Gene Expression Before and After Ambient RNA Correction"
MARKER_TSNE_GENE_TITLE_SIZE = 14
MARKER_TSNE_LABEL_SIZE = 13
MARKER_TSNE_MAIN_TITLE_SIZE = 17
MARKER_TSNE_MAX_PER_ROW = 4
MARKER_TSNE_BEFORE_AFTER_HSPACE = 0.26
MARKER_TSNE_BAND_HSPACE = 0.38
MARKER_TSNE_COLUMN_WSPACE = 0.50
MARKER_TSNE_DPI = 150


def _marker_gene_bands(genes: list[str]) -> list[list[str]]:
    """Chunk genes into bands of at most MARKER_TSNE_MAX_PER_ROW, stacked vertically."""
    return [
        genes[i: i + MARKER_TSNE_MAX_PER_ROW]
        for i in range(0, len(genes), MARKER_TSNE_MAX_PER_ROW)
    ]


def _draw_marker_tsne_cell(
    ax,
    coords: np.ndarray,
    vals: np.ndarray,
    *,
    gene: str,
    row_label: str | None,
) -> None:
    import matplotlib.pyplot as plt

    sc = ax.scatter(
        coords[:, 0], coords[:, 1], c=vals, cmap="Purples",
        s=6, alpha=0.8, linewidths=0, rasterized=True,
    )
    plt.colorbar(sc, ax=ax, pad=0.04, fraction=0.040)
    ax.set_title(gene, fontsize=MARKER_TSNE_GENE_TITLE_SIZE, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if row_label is not None:
        ax.set_ylabel(
            row_label,
            fontsize=MARKER_TSNE_LABEL_SIZE,
            fontweight="bold",
        )


def plot_marker_genes_tsne(
    tsne_before: np.ndarray,
    tsne_after: np.ndarray,
    expr_before: "dict[str, np.ndarray]",
    expr_after: "dict[str, np.ndarray]",
    *,
    out_dir: "Path | str",
    stem: str,
    genes: list[str] | None = None,
) -> "list[Path]":
    """Before/after marker expression on t-SNE; stacks bands of up to 4 genes."""
    _apply_style()
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_theme(style="ticks")

    if genes is None:
        plot_genes = [g for g in expr_before if g in expr_after]
    else:
        plot_genes = [g for g in genes if g in expr_before and g in expr_after]
    if not plot_genes:
        return []

    bands = _marker_gene_bands(plot_genes)
    n_bands = len(bands)
    fig_w = 4.5 * MARKER_TSNE_MAX_PER_ROW
    fig_h = 4.2 * n_bands + 1.0
    fig = plt.figure(figsize=(fig_w, fig_h))
    outer = gridspec.GridSpec(n_bands, 1, figure=fig, hspace=MARKER_TSNE_BAND_HSPACE)

    tsne_before = np.asarray(tsne_before)
    tsne_after = np.asarray(tsne_after)
    row_specs = (
        (tsne_before, expr_before, "Before"),
        (tsne_after, expr_after, "After"),
    )

    for band_idx, band_genes in enumerate(bands):
        inner = gridspec.GridSpecFromSubplotSpec(
            2, MARKER_TSNE_MAX_PER_ROW,
            subplot_spec=outer[band_idx],
            hspace=MARKER_TSNE_BEFORE_AFTER_HSPACE,
            wspace=MARKER_TSNE_COLUMN_WSPACE,
        )
        for col, gene in enumerate(band_genes):
            for row, (coords, expr_map, row_label) in enumerate(row_specs):
                ax = fig.add_subplot(inner[row, col])
                vals = np.asarray(expr_map[gene], dtype=float)
                _draw_marker_tsne_cell(
                    ax, coords, vals,
                    gene=gene,
                    row_label=row_label if col == 0 else None,
                )

    fig.suptitle(
        MARKER_TSNE_MAIN_TITLE,
        fontsize=MARKER_TSNE_MAIN_TITLE_SIZE,
        fontweight="bold",
        y=0.995,
    )
    fig.subplots_adjust(top=0.945)
    # PNG only at moderate DPI — large multi-panel figures OOM on login nodes when
    # saving 300 dpi PNG + vector PDF together (22+ subplots for 11 genes).
    return save_figure(
        fig, out_dir, stem, also_pdf=False, dpi=MARKER_TSNE_DPI,
    )


def plot_counts_before_after(pre: np.ndarray, post: np.ndarray, *,
                              out_dir: Path | str, stem: str, title: str) -> list[Path]:
    """Scatter of per-cell total counts before vs after ambient correction."""
    _apply_style()
    import matplotlib.pyplot as plt
    pre = np.asarray(pre, dtype=float)
    post = np.asarray(post, dtype=float)
    finite = np.isfinite(pre) & np.isfinite(post) & (pre > 0)
    pre, post = pre[finite], post[finite]
    fig, ax = plt.subplots(figsize=(6.0, 5.5))
    if pre.size == 0:
        ax.set_title(title + " (no data)")
        return save_figure(fig, out_dir, stem)
    ax.scatter(pre, post, s=6, alpha=0.4, color="#1f77b4", linewidths=0)
    lim_max = float(max(pre.max(), post.max()))
    ax.plot([0, lim_max], [0, lim_max], color="black", linewidth=0.8,
            linestyle="--", label="y = x (no correction)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(left=max(1.0, float(pre.min())))
    ax.set_ylim(bottom=max(1.0, float(post[post > 0].min()) if (post > 0).any() else 1.0))
    ax.set_xlabel("total counts (pre-correction)")
    ax.set_ylabel("total counts (post-correction)")
    ax.set_title(title)
    ax.legend(fontsize=FONT_SIZE - 1, loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return save_figure(fig, out_dir, stem)
