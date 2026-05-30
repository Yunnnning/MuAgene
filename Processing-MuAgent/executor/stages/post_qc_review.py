"""post_qc_review — QC review user checkpoint (checkpoint #2 of 3).

Runs after S3 doublet filtering and before S4/S5 dimensionality reduction.
Generates QC figures and writes the QC review summary at
  deliverables/checkpoint/qc_review/qc_summary.md

This is the single user-facing QC checkpoint: inspect S1/S2 QC figures, S3
doublet histograms, and the cell-count waterfall; adjust S1/S2 thresholds or
(paired only) the S3 union doublet policy if needed; then approve to
continue. S3 executes before this stage — paired doublet policy is surfaced
here, not at a separate S3 gate.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .. import qc_summary as _qcs
from ..log import log_event
from ..run_paths import RunPaths


def _plot_score_hist(
    scores: np.ndarray,
    flags: np.ndarray,
    *,
    title: str,
    out_dir: Path,
    stem: str,
) -> list[Path]:
    from .. import figures as _fig
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    finite = np.isfinite(scores)
    scores, flags = scores[finite], flags[finite]
    if scores.size == 0:
        ax.set_title(title + " (no data)")
        return _fig.save_figure(fig, out_dir, stem)

    bins = np.linspace(0.0, max(1.0, float(scores.max()) * 1.05), 51)
    n_sing = int((~flags).sum())
    n_doub = int(flags.sum())
    ax.hist(scores[~flags], bins=bins, color="#3b82f6", alpha=0.75,
            label=f"singlet (n={n_sing})", edgecolor="none")
    ax.hist(scores[flags], bins=bins, color="#ef4444", alpha=0.75,
            label=f"doublet (n={n_doub})", edgecolor="none")
    ax.set_xlabel("doublet score")
    ax.set_ylabel("cells")
    ax.set_title(title)
    ax.legend(fontsize=_fig.FONT_SIZE - 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return _fig.save_figure(fig, out_dir, stem)


def _plot_doublet_scores(run_dir: Path, figs_dir: Path) -> list[Path]:
    """Doublet-score histograms (RNA + ATAC) from S3 calls.parquet."""
    calls_path = RunPaths(run_dir).stage_dir("s3_doublets") / "calls.parquet"
    if not calls_path.exists():
        return []

    calls = pd.read_parquet(calls_path)
    result: list[Path] = []

    if "scrublet_score" in calls.columns:
        rna_scores = calls["scrublet_score"].dropna().to_numpy().astype(float)
        rna_flags = (
            calls.loc[calls["scrublet_score"].notna(), "scrublet_is_doublet"]
            .fillna(False).to_numpy().astype(bool)
        )
        if rna_scores.size > 0:
            result.extend(_plot_score_hist(
                rna_scores, rna_flags,
                title="RNA doublet scores (Scrublet)",
                out_dir=figs_dir,
                stem="post_qc_review_doublet_rna",
            ))

    if "atac_doublet_score" in calls.columns:
        atac_scores = calls["atac_doublet_score"].dropna().to_numpy().astype(float)
        atac_flags = (
            calls.loc[calls["atac_doublet_score"].notna(), "atac_is_doublet"]
            .fillna(False).to_numpy().astype(bool)
        )
        if atac_scores.size > 0:
            result.extend(_plot_score_hist(
                atac_scores, atac_flags,
                title="ATAC doublet scores (SnapATAC2)",
                out_dir=figs_dir,
                stem="post_qc_review_doublet_atac",
            ))

    return result


def _plot_cell_count_waterfall(run_dir: Path, figs_dir: Path) -> list[Path]:
    """Grouped bar chart: RNA and ATAC cell counts at each preprocessing stage."""
    from .. import figures as _fig
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    counts = _qcs._stage_counts(run_dir)

    rna_ingest = counts.get("rna_after_ambient") or counts.get("rna_ingest")
    stages = [
        ("raw",          counts.get("rna_raw"),          counts.get("atac_raw_barcodes")),
        ("after ambient\n/ ingest", rna_ingest,          counts.get("atac_raw_barcodes")),
        ("after\nS1/S2 QC",  counts.get("rna_qc_post"),  counts.get("atac_qc_post")),
        ("after\nS3 doublets", counts.get("rna_post_doublet"), counts.get("atac_post_doublet")),
    ]

    stages = [(lbl, r, a) for lbl, r, a in stages if r is not None or a is not None]
    if not stages:
        return []

    labels = [s[0] for s in stages]
    rna_h = [s[1] if s[1] is not None else 0 for s in stages]
    atac_h = [s[2] if s[2] is not None else 0 for s in stages]

    x = np.arange(len(labels))
    width = 0.35

    _fig._apply_style()
    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    bars_rna = ax.bar(x - width / 2, rna_h, width, label="RNA", color="#3b82f6", alpha=0.85)
    bars_atac = ax.bar(x + width / 2, atac_h, width, label="ATAC", color="#f97316", alpha=0.85)
    ax.bar_label(
        bars_rna,
        labels=[f"{v:,}" if v > 0 else "" for v in rna_h],
        padding=3,
        fontsize=_fig.FONT_SIZE - 2,
    )
    ax.bar_label(
        bars_atac,
        labels=[f"{v:,}" if v > 0 else "" for v in atac_h],
        padding=3,
        fontsize=_fig.FONT_SIZE - 2,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, ha="center")
    ax.set_ylabel("cells")
    ax.set_title("Cell count across preprocessing stages (S0–S3)")
    ax.legend()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return _fig.save_figure(fig, figs_dir, "post_qc_review_cell_counts")


def propose(run_dir: Path | str) -> dict[str, Any]:
    """Generate QC figures and summary markdown; return proposal content dict."""
    run_dir = Path(run_dir)
    rp = RunPaths(run_dir)
    figs_dir = rp.deliv_qc_review
    figs_dir.mkdir(parents=True, exist_ok=True)

    figures_generated: list[str] = []

    try:
        paths = _plot_cell_count_waterfall(run_dir, figs_dir)
        figures_generated.extend(str(p) for p in paths)
    except Exception as e:
        log_event(run_dir, {"stage": "post_qc_review", "event": "waterfall_failed",
                             "error": str(e)})

    try:
        paths = _plot_doublet_scores(run_dir, figs_dir)
        figures_generated.extend(str(p) for p in paths)
    except Exception as e:
        log_event(run_dir, {"stage": "post_qc_review", "event": "doublet_hist_failed",
                             "error": str(e)})

    summary_path = rp.qc_review_summary_md
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        summary_path.write_text(_qcs.build_qc_review(run_dir))
    except Exception as e:
        log_event(run_dir, {"stage": "post_qc_review", "event": "summary_failed",
                             "error": str(e)})
        summary_path.write_text(f"# QC review checkpoint\n\n_Error building summary: {e}_\n")

    log_event(run_dir, {
        "stage": "post_qc_review",
        "event": "propose_done",
        "n_figures": len(figures_generated),
        "summary": str(summary_path),
    })
    return {
        "figures": figures_generated,
        "qc_summary": str(summary_path),
    }
