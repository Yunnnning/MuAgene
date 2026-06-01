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
TITLE_SIZE = 14
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


def plot_fragment_size_distribution(distr: np.ndarray, *, out_dir: Path | str,
                                     stem: str, title: str) -> list[Path]:
    """Plot a 1D fragment-size histogram (counts per fragment length).

    `distr` is a 1D vector where `distr[i]` is the number of fragments of length
    `i` (the SnapATAC2 `frag_size_distr` format; index 0 holds the over-max bucket).
    A well-prepared ATAC library shows distinct peaks at ~150, ~300, ~450 bp.
    """
    _apply_style()
    import matplotlib.pyplot as plt
    d = np.asarray(distr, dtype=float).ravel()
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    if d.size == 0:
        ax.set_title(title + " (no data)")
        return save_figure(fig, out_dir, stem)
    # Drop the over-max bucket (index 0) — its size is dataset-dependent and
    # would dominate the y-axis on libraries with many large fragments.
    body = d[1:]
    x = np.arange(1, body.size + 1)
    ax.fill_between(x, 0, body, alpha=0.55, color="#1f77b4", linewidth=0)
    ax.plot(x, body, color="#1f4f8b", linewidth=1.0)
    x_right = min(1000, body.size)
    ax.set_xlim(left=0, right=x_right)
    y_top = ax.get_ylim()[1]
    for vline in (147, 294, 441):
        if vline < body.size:
            ax.axvline(vline, color="firebrick", linestyle=":", linewidth=0.9, alpha=0.6)
    for x0, x1, label in [(1, 147, "nucleosome-free"),
                          (147, 294, "mono"),
                          (294, 441, "di"),
                          (441, x_right, "tri")]:
        xc = (x0 + x1) / 2
        if xc < x_right:
            ax.text(xc, y_top * 0.97, label, ha="center", va="top",
                    fontsize=FONT_SIZE - 3, color="firebrick")
    ax.set_xlabel("fragment length (bp)")
    ax.set_ylabel("fragment count")
    ax.set_title(title)
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
