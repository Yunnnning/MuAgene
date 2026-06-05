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
FONT_SIZE = 12
MEDIUM_BLUE = "mediumblue"
QC_FILL_COLOR = "#f97316"
QC_FILL_ALPHA = 0.50
QC_EDGE_COLOR = "#c2410c"
QC_EDGE_LINEWIDTH = 0.4
ANNOTATION_LINEWIDTH = 1.5
FSD_ANNOTATION_LINEWIDTH = 1.2
ANNOTATION_FONTSIZE = FONT_SIZE - 2
TITLE_SIZE = 14
FRIP_DISTRIBUTION_TITLE = "Fraction of Reads in Peaks (FRiP) distribution"
UMAP_SIZE = (6.5, 5.5)
QC_VIOLIN_SIZE = (12, 4.5)


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


def save_figure(fig, out_dir: Path | str, stem: str, *, also_pdf: bool = True) -> list[Path]:
    """Save `fig` as <stem>.png (300 dpi) and optionally <stem>.pdf. Returns list of paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    png_path = out_dir / f"{stem}.png"
    fig.savefig(png_path, dpi=FIGURE_DPI, bbox_inches="tight")
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
        if v.size == 0:
            ax.set_title(name + " (no data)")
            continue
        parts = ax.violinplot(v, showmeans=False, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_alpha(0.65)
        ax.set_title(name)
        ax.set_ylabel(name)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(title)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
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

    def _prep(d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        d = np.asarray(d, dtype=float).ravel()
        body = d[1:]  # drop over-max bucket at index 0
        total = body.sum()
        normed = body / total if total > 0 else body
        x = np.arange(1, body.size + 1)
        return x, normed

    d = np.asarray(distr, dtype=float).ravel()
    if d.size == 0:
        fig, ax = plt.subplots(figsize=(8.0, 4.5))
        ax.set_title(title + " (no data)")
        return save_figure(fig, out_dir, stem)

    x_b, body_b = _prep(d)
    x_right = min(1000, body_b.size)
    ylabel = "fraction of fragments"

    def _annotate_nucleosome(ax, y_top: float) -> None:
        for vline in (147, 294, 441):
            if vline < x_right:
                ax.axvline(
                    vline, color=MEDIUM_BLUE, linestyle="--",
                    linewidth=FSD_ANNOTATION_LINEWIDTH, zorder=5,
                )
        for x0, x1, label in [(1, 147, "nucleosome\nfree"),
                              (147, 294, "mono"),
                              (294, 441, "di"),
                              (441, x_right, "tri")]:
            xc = (x0 + x1) / 2
            if xc < x_right:
                ax.text(
                    xc, y_top * 0.97, label, ha="center", va="top",
                    fontsize=ANNOTATION_FONTSIZE, color=MEDIUM_BLUE,
                )

    if distr_after is not None and np.asarray(distr_after).size > 0:
        x_plot, body_plot = _prep(np.asarray(distr_after, dtype=float))
    else:
        x_plot, body_plot = x_b, body_b

    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    x = x_plot[:x_right]
    y = body_plot[:x_right]
    ax.fill_between(x, 0, y, alpha=QC_FILL_ALPHA, color=QC_FILL_COLOR, linewidth=0)
    ax.plot(x, y, color=QC_EDGE_COLOR, linewidth=QC_EDGE_LINEWIDTH)
    ax.set_xlim(left=0, right=x_right)
    ax.set_ylabel(ylabel)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    y_top = float(y.max()) if y.size else 1.0
    _annotate_nucleosome(ax, y_top)
    ax.set_xlabel("fragment length (bp)")
    ax.set_title(title)

    return save_figure(fig, out_dir, stem)


def plot_frip_histogram(
    frip: np.ndarray,
    *,
    out_dir: Path | str,
    stem: str,
    frip_min: float = 0.2,
) -> list[Path]:
    """Histogram of per-cell FRiP values with a vertical threshold line.

    `frip` contains FRiP values for all cells that passed the 3-metric QC
    (before the FRiP filter itself). The vertical line shows the chosen threshold.
    """
    _apply_style()
    import matplotlib.pyplot as plt

    frip = np.asarray(frip, dtype=float)
    frip = frip[np.isfinite(frip)]

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    if frip.size == 0:
        ax.set_title(f"{FRIP_DISTRIBUTION_TITLE} (no data)")
        return save_figure(fig, out_dir, stem)

    n_total = frip.size
    n_pass = int((frip >= frip_min).sum())
    pct_pass = 100.0 * n_pass / n_total if n_total > 0 else 0.0

    bins = np.linspace(0, 1, 51)
    ax.hist(
        frip, bins=bins, color=QC_FILL_COLOR, alpha=QC_FILL_ALPHA,
        edgecolor=QC_EDGE_COLOR, linewidth=QC_EDGE_LINEWIDTH,
    )
    ax.axvline(
        frip_min, color=MEDIUM_BLUE, linestyle="--",
        linewidth=ANNOTATION_LINEWIDTH, zorder=5,
    )
    ax.text(
        0.98, 0.98,
        f"threshold = {frip_min:.2f}\n{n_pass}/{n_total} passed ({pct_pass:.1f}%)",
        transform=ax.transAxes,
        ha="right", va="top", fontsize=ANNOTATION_FONTSIZE, color=MEDIUM_BLUE,
    )
    ax.set_xlabel("FRiP (fraction of reads in peaks)")
    ax.set_ylabel("number of cells")
    ax.set_title(FRIP_DISTRIBUTION_TITLE)
    ax.set_xlim(0, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return save_figure(fig, out_dir, stem)


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
