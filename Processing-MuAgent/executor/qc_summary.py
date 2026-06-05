"""Concise user-facing QC summary (markdown).

Audit-driven rewrite — reporting-only changes, no pipeline-logic edits:
    1. Cell-count flow table across stages (makes all transitions visible).
    2. Hidden SnapATAC2 import-stage cell drop is surfaced explicitly.
    3. Doublet overlap table restricted to cells evaluated by BOTH detectors;
       "not evaluated" is reported separately and no longer conflated with
       "not flagged".
    4. Per-modality doublet removal counts (computed from post-doublet h5ads).
    5. Paired-intersection happens at S3 (not S8); n_cells_joint is surfaced
       in the flow table and "Final retained" section. S8 assembly is a safety
       no-op intersection.
    6. Baselines relabelled from ambiguous "Cells before filtering" to
       stage-aware phrasing (e.g. "Cells entering this stage").
    7. Thresholds rounded to 2 decimal places; integer-typed values are
       displayed as integers.
"""
from __future__ import annotations

import base64
import html as html_module
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


@dataclass(frozen=True)
class _FigRender:
    """Figure inclusion and relative paths for markdown/HTML output."""

    md_parent: Path | None = None
    embed_figures: bool = True  # True for qc_review.md and HTML report; False for raw-text excerpts


_DEFAULT_FIG_RENDER = _FigRender()


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _pct(a: int, b: int) -> str:
    return f"{(100.0 * a / b):.1f}%" if b else "n/a"


def _fmt(value: Any) -> str:
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
        return ", ".join(_fmt(x) for x in value)
    return str(value)


def _md_table_cell(value: Any) -> str:
    """Format a table cell; preserve inline markdown (**, trailing *)."""
    if isinstance(value, str) and ("**" in value or value.endswith("*")):
        return value
    return _fmt(value)


def _md_table(header: list[str], rows: list[list[Any]]) -> str:
    align = "|" + "|".join("---" for _ in header) + "|"
    h = "| " + " | ".join(_md_table_cell(c) for c in header) + " |"
    body = "\n".join(
        "| " + " | ".join(_md_table_cell(x) for x in r) + " |" for r in rows
    )
    return f"{h}\n{align}\n{body}"


def _stats_row(name: str, vals: np.ndarray) -> list[Any]:
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return [name, "n/a", "n/a", "n/a", "n/a"]
    return [
        name,
        f"{np.mean(vals):.2f}",
        f"{np.median(vals):.2f}",
        f"{np.min(vals):.2f}",
        f"{np.max(vals):.2f}",
    ]


def _param(params: dict[str, Any], key: str) -> Any:
    entry = params.get(key)
    return entry.get("value") if isinstance(entry, dict) else None


def _fig_block(
    run_dir: Path,
    stem: str,
    *,
    caption: str | None = None,
    render: _FigRender = _DEFAULT_FIG_RENDER,
) -> str:
    """Embed a checkpoint QC figure (markdown/HTML) or skip when embed_figures=False."""
    from .run_paths import RunPaths
    if not render.embed_figures:
        return ""
    png = RunPaths(run_dir).deliv_qc_review / f"{stem}.png"
    if not png.exists():
        return ""
    alt = caption or stem.replace("_", " ")
    if render.md_parent is not None:
        src = Path(os.path.relpath(png, render.md_parent)).as_posix()
    else:
        src = f"./{stem}.png"
    return f"\n![{alt}]({src})\n"


def _fig_pair_block(
    run_dir: Path,
    left: tuple[str, str],
    right: tuple[str, str],
    *,
    render: _FigRender = _DEFAULT_FIG_RENDER,
) -> str:
    """Embed two checkpoint QC figures side-by-side (HTML in markdown reports)."""
    from .run_paths import RunPaths
    if not render.embed_figures:
        return ""
    fig_dir = RunPaths(run_dir).deliv_qc_review
    panels: list[str] = []
    for stem, caption in (left, right):
        png = fig_dir / f"{stem}.png"
        if not png.exists():
            continue
        if render.md_parent is not None:
            src = Path(os.path.relpath(png, render.md_parent)).as_posix()
        else:
            src = f"./{stem}.png"
        alt = html_module.escape(caption)
        panels.append(
            '<figure class="qc-figure" style="flex:1 1 45%; min-width:260px; margin:0; '
            'display:flex; flex-direction:column;">'
            f'<img src="{src}" alt="{alt}" style="width:100%; height:300px; '
            'object-fit:contain; object-position:center top; display:block;">'
            f'<figcaption style="font-size:0.85rem; color:#555; margin-top:0.35rem;">'
            f"{alt}</figcaption></figure>"
        )
    if len(panels) < 2:
        return ""
    return (
        '\n<div class="qc-plot-pair" style="display:flex; flex-wrap:wrap; gap:1rem; '
        'align-items:stretch; margin:1rem 0;">'
        + "".join(panels)
        + "</div>\n"
    )


def _inline_md(text: str) -> str:
    """Bold/italic on already-escaped HTML text (preserve snake_case names)."""
    s = html_module.escape(text)
    s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"<em>\1</em>", s)
    return s


def _md_table_html(table_lines: list[str]) -> str:
    rows_html: list[str] = []
    for row in table_lines:
        if re.fullmatch(r"\|[\s\-:|]+\|", row.strip()):
            continue
        cells = [c.strip() for c in row.strip().strip("|").split("|")]
        if not rows_html:
            rows_html.append(
                "<tr>" + "".join(f"<th>{_inline_md(c)}</th>" for c in cells) + "</tr>"
            )
        else:
            rows_html.append(
                "<tr>" + "".join(f"<td>{_inline_md(c)}</td>" for c in cells) + "</tr>"
            )
    return "<table>" + "".join(rows_html) + "</table>"


def _markdown_to_html(text: str) -> str:
    """Convert QC-summary markdown subset to HTML (no third-party deps)."""
    lines = text.splitlines()
    out: list[str] = []
    i = 0
    img_re = re.compile(r"^!\[([^\]]*)\]\(([^)]+)\)\s*$")
    hdr_re = re.compile(r"^(#{1,6})\s+(.*)$")

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if stripped.startswith("|"):
            block: list[str] = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                block.append(lines[i])
                i += 1
            out.append(_md_table_html(block))
            continue

        hdr = hdr_re.match(line)
        if hdr:
            level = len(hdr.group(1))
            out.append(f"<h{level}>{_inline_md(hdr.group(2))}</h{level}>")
            i += 1
            continue

        img = img_re.match(stripped)
        if img:
            alt, src = img.group(1), img.group(2)
            out.append(
                f'<p><img src="{html_module.escape(src, quote=True)}" '
                f'alt="{html_module.escape(alt)}"></p>'
            )
            i += 1
            continue

        if not stripped:
            i += 1
            continue

        if stripped.startswith("- "):
            items: list[str] = []
            while i < len(lines) and lines[i].strip().startswith("- "):
                item = lines[i].strip()[2:]
                cls = ' class="qc-footnote"' if item.startswith("*") else ""
                items.append(f"<li{cls}>{_inline_md(item)}</li>")
                i += 1
            out.append("<ul>" + "".join(items) + "</ul>")
            continue

        # Footnote annotation (line begins with * referencing table markers like n=7200*)
        if stripped.startswith("*") and not stripped.startswith("**"):
            out.append(f'<p class="qc-footnote">{_inline_md(stripped)}</p>')
            i += 1
            continue

        para = [_inline_md(line)]
        i += 1
        while i < len(lines):
            nxt = lines[i].strip()
            if (
                not nxt
                or nxt.startswith("|")
                or hdr_re.match(lines[i])
                or nxt.startswith("- ")
                or img_re.match(nxt)
            ):
                break
            para.append(_inline_md(lines[i]))
            i += 1
        out.append("<p>" + "<br>\n".join(para) + "</p>")

    return "\n".join(out)


def _embed_html_images(html: str, fig_dir: Path) -> str:
    """Inline sibling PNGs as data URIs so the report works when opened anywhere."""

    def _repl(match: re.Match[str]) -> str:
        src, alt = match.group(1), match.group(2)
        if src.startswith("data:"):
            return match.group(0)
        path = fig_dir / src
        if not path.is_file():
            return match.group(0)
        data = base64.standard_b64encode(path.read_bytes()).decode("ascii")
        return (
            f'<img src="data:image/png;base64,{data}" '
            f'alt="{html_module.escape(alt)}">'
        )

    return re.sub(
        r'<img src="([^"]+)" alt="([^"]*)">',
        _repl,
        html,
    )


# ---------------------------------------------------------------------------
# Stage probes — read current state of each stage output
# ---------------------------------------------------------------------------

def _stage_counts(run_dir: Path) -> dict[str, Any]:
    """Collect cell counts at every meaningful transition. Returns dict with
    keys for each stage; missing stages get None."""
    from .run_paths import RunPaths
    A = RunPaths(run_dir).artifacts
    counts: dict[str, Any] = {
        "rna_raw": None, "atac_raw_barcodes": None,
        "rna_ingest": None,
        "rna_after_ambient": None,
        "atac_after_snap_import": None,
        "rna_qc_post": None, "atac_qc_post": None,
        "rna_post_doublet": None, "atac_post_doublet": None,
        "n_cells_joint": None,
        "rna_final": None, "atac_final": None,
    }

    # S0 validation report
    vr_path = A / "s0_ingest" / "validation_report.json"
    if vr_path.exists():
        vr = json.loads(vr_path.read_text())
        counts["rna_raw"] = int(vr.get("rna_n_cells", 0))
        counts["atac_raw_barcodes"] = int(vr.get("atac_n_unique_barcodes", 0))

    # S1a ambient correction output
    rna_ambient = A / "s1a_ambient" / "rna_decontaminated.h5ad"
    if rna_ambient.exists():
        try:
            import anndata as ad
            a = ad.read_h5ad(rna_ambient, backed="r")
            counts["rna_after_ambient"] = int(a.n_obs)
            try: a.file.close()
            except Exception: pass
        except Exception:
            pass

    # S0 ingest paired RNA
    rna_ingest = A / "s0_ingest" / "rna_ingest.h5ad"
    if rna_ingest.exists():
        try:
            import anndata as ad
            a = ad.read_h5ad(rna_ingest, backed="r")
            counts["rna_ingest"] = int(a.n_obs)
            try: a.file.close()
            except Exception: pass
        except Exception:
            pass

    # S2 atac_qc summary reports counts around the SnapATAC2 import
    atac_summary = A / "s2_atac_qc" / "qc_summary.json"
    if atac_summary.exists():
        s2 = json.loads(atac_summary.read_text())
        counts["atac_after_snap_import"] = int(s2.get("n_cells_pre", 0))
        counts["atac_qc_post"] = int(s2.get("n_cells_post", 0))

    # S1 RNA QC post
    s1_post = A / "s1_rna_qc" / "qc_metrics_post.parquet"
    if s1_post.exists():
        counts["rna_qc_post"] = int(len(pd.read_parquet(s1_post)))

    # S3 doublet outputs
    rna_pd = A / "s3_doublets" / "rna_post_doublet.h5ad"
    if rna_pd.exists():
        try:
            import anndata as ad
            a = ad.read_h5ad(rna_pd, backed="r")
            counts["rna_post_doublet"] = int(a.n_obs)
            try: a.file.close()
            except Exception: pass
        except Exception:
            pass
    atac_pd = A / "s3_doublets" / "atac_post_doublet.h5ad"
    if atac_pd.exists():
        try:
            import snapatac2 as snap
            a = snap.read(str(atac_pd))
            counts["atac_post_doublet"] = int(a.n_obs)
            try: a.close()
            except Exception: pass
        except Exception:
            pass

    # S3 paired-intersection sentinel — present on the paired branch only.
    joint_path = A / "s3_doublets" / "joint_barcodes.txt"
    if joint_path.exists():
        try:
            text = joint_path.read_text()
            counts["n_cells_joint"] = sum(1 for line in text.splitlines() if line.strip())
        except Exception:
            pass

    # S8 final
    processed = A / "s8_umap" / "processed.h5mu"
    rna_h5ad_final = A / "s8_umap" / "rna_processed.h5ad"
    atac_h5ad_final = A / "s8_umap" / "atac_processed.h5ad"
    try:
        import mudata as mu
        import anndata as ad
        if processed.exists():
            m = mu.read_h5mu(str(processed))
            counts["rna_final"] = int(m.mod["rna"].n_obs) if "rna" in m.mod else None
            counts["atac_final"] = int(m.mod["atac"].n_obs) if "atac" in m.mod else None
        elif rna_h5ad_final.exists():
            counts["rna_final"] = int(ad.read_h5ad(rna_h5ad_final, backed="r").n_obs)
            counts["atac_final"] = int(ad.read_h5ad(atac_h5ad_final, backed="r").n_obs)
    except Exception:
        pass

    return counts


def _barcode_set_h5ad(path: Path) -> set[str] | None:
    if not path.exists():
        return None
    try:
        import anndata as ad
        a = ad.read_h5ad(path, backed="r")
        bc = set(a.obs_names)
        try:
            a.file.close()
        except Exception:
            pass
        return bc
    except Exception:
        return None


def _barcode_set_snap(path: Path) -> set[str] | None:
    if not path.exists():
        return None
    try:
        import snapatac2 as snap
        a = snap.read(str(path))
        bc = set(a.obs_names)
        try:
            a.close()
        except Exception:
            pass
        return bc
    except Exception:
        return None


def _paired_shared_flow_counts(run_dir: Path, counts: dict[str, Any]) -> list[int | None]:
    """Shared RNA∩ATAC barcode counts per flow row (paired branch only)."""
    from .run_paths import RunPaths
    A = RunPaths(run_dir).artifacts
    n_rows = 5 if counts.get("rna_final") is not None else 4
    shared: list[int | None] = [None] * n_rows

    rna_ingest = _barcode_set_h5ad(A / "s0_ingest" / "rna_ingest.h5ad")
    if rna_ingest is None:
        return shared

    ingest_shared = len(rna_ingest)
    shared[0] = ingest_shared
    shared[1] = ingest_shared

    calls_p = A / "s3_doublets" / "calls.parquet"
    if calls_p.exists():
        calls = pd.read_parquet(calls_p)
        both = calls["scrublet_score"].notna() & calls["atac_doublet_score"].notna()
        shared[2] = int(both.sum())

    joint = counts.get("n_cells_joint")
    if joint is not None:
        shared[3] = int(joint)
    if n_rows > 4 and counts.get("rna_final") is not None:
        shared[4] = int(counts["rna_final"])

    return shared


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------

def _flow_section(
    run_dir: Path,
    counts: dict[str, Any],
    workflow_branch: str,
    *,
    include_final_stage: bool = True,
) -> str:
    """Cell-count flow across stages; makes transitions visible end-to-end.

    QC review checkpoint summaries stop at S3 (``include_final_stage=False``);
    the post-run manifest summary includes S4–S8 when those stages have run.
    """
    def fmt(v):
        return "n/a" if v is None else str(int(v))

    rna_raw = counts["rna_raw"]
    atac_raw = counts["atac_raw_barcodes"]

    joint = counts.get("n_cells_joint")

    rna_after_s1a = counts["rna_after_ambient"]
    if rna_after_s1a is None:
        rna_after_s1a = counts["rna_ingest"]

    s3_note = (
        "union doublet removal and joint cell retention"
        if joint is not None
        else "union doublet removal per modality"
    )

    paired = workflow_branch == "paired"
    shared_counts = _paired_shared_flow_counts(run_dir, counts) if paired else []

    def _flow_row(stage: str, rna: Any, atac: Any, note: str, shared_idx: int) -> list[Any]:
        row = [stage, fmt(rna), fmt(atac)]
        if paired:
            sh = shared_counts[shared_idx] if shared_idx < len(shared_counts) else None
            row.append(fmt(sh))
        row.append(note)
        return row

    rows = [
        _flow_row("1. raw", rna_raw, atac_raw, "—", 0),
        _flow_row("2. after ambient RNA correction", rna_after_s1a, atac_raw,
                  "RNA ambient correction (cell count unchanged)", 1),
        _flow_row("3. after RNA / ATAC QC", counts["rna_qc_post"], counts["atac_qc_post"],
                  "per-modality MAD and quality thresholds", 2),
        _flow_row("4. after doublet removal", counts["rna_post_doublet"], counts["atac_post_doublet"],
                  s3_note, 3),
    ]
    if include_final_stage:
        rna_inter = (
            (counts["rna_post_doublet"] - counts["rna_final"])
            if (counts["rna_post_doublet"] is not None and counts["rna_final"] is not None) else None
        )
        atac_inter = (
            (counts["atac_post_doublet"] - counts["atac_final"])
            if (counts["atac_post_doublet"] is not None and counts["atac_final"] is not None) else None
        )
        rows.append(
            _flow_row("5. after S4–S8 (final)", counts["rna_final"], counts["atac_final"],
                      (f"paired: S8 assembly is a no-op intersection; RNA lost {fmt(rna_inter)}, "
                       f"ATAC lost {fmt(atac_inter)} downstream of S3"
                       if joint is not None
                       else "per-modality final outputs (no joint object on this branch)"),
                      4)
        )
    headers = ["stage", "RNA", "ATAC"] + (["Shared"] if paired else []) + ["note"]
    return (
        "## Cell-count flow across stages\n\n"
        f"{_md_table(headers, rows)}\n"
    )


def _flow_figures(run_dir: Path, *, render: _FigRender = _DEFAULT_FIG_RENDER) -> str:
    return _fig_block(
        run_dir, "post_qc_review_cell_counts",
        caption="Cell counts across preprocessing (RNA and ATAC).",
        render=render,
    )


def _ambient_section(
    run_dir: Path,
    params: dict[str, Any],
    counts: dict[str, Any],
    *,
    render: _FigRender = _DEFAULT_FIG_RENDER,
) -> str:
    from .run_paths import RunPaths
    s1a = RunPaths(run_dir).stage_dir("s1a_ambient")
    summary_p = s1a / "summary.json"

    method = _param(params, "s1a_ambient.method")
    if method in (None, "none", "skipped_empty", "skipped_no_r"):
        note = {
            None: "_(stage did not run; legacy run or RNA absent)_",
            "none": "_disabled in preprocessing plan (method=none); "
                    "no ambient correction applied._",
            "skipped_empty": "_RNA AnnData empty (atac_only branch); pass-through._",
            "skipped_no_r": (
                "**SKIPPED — R packages unavailable at runtime.** "
                "Counts are uncorrected (pass-through). "
                "Install missing R packages and re-run ambient RNA correction before approving QC."
            ),
        }[method]
        return "## Ambient RNA correction\n\n" + note + "\n"

    summary: dict[str, Any] = {}
    if summary_p.exists():
        try:
            summary = json.loads(summary_p.read_text())
        except Exception:
            summary = {}

    median_rho = (
        summary.get("median_contamination")
        if summary.get("median_contamination") is not None
        else _param(params, "s1a_ambient.median_contamination")
    )
    pre_total = summary.get("total_counts_pre")
    post_total = summary.get("total_counts_post")
    pct_removed = ""
    if pre_total and post_total and int(pre_total) > 0:
        pct_removed = f" ({_pct(int(pre_total) - int(post_total), int(pre_total))} of UMIs removed)"

    rho_note = ""
    if median_rho is not None:
        rho_note = (
            f"**rho** (median {_fmt(median_rho)}): estimated fraction of each cell's "
            "counts attributed to ambient RNA before correction "
            "(SoupX/DecontX often apply one global or cluster-level estimate to all cells).\n"
        )

    counts_note = ""
    if pre_total is not None and post_total is not None:
        counts_note = (
            f"**Total UMI counts** (sum over all cells): {_fmt(pre_total)} pre-correction → "
            f"{_fmt(post_total)} post-correction{pct_removed}.\n"
        )

    rows = [
        ["method", method],
        ["rho (median)", median_rho],
        ["high-contamination cells (rho>0.20)",
         _param(params, "s1a_ambient.n_high_contamination_cells")],
        ["max-contamination cap",
         summary.get("max_contam_cap", _param(params, "s1a_ambient.max_contamination"))],
        ["pre-correction total counts (sum)", pre_total],
        ["post-correction total counts (sum)", post_total],
    ]

    fig = _fig_block(
        run_dir, "s1a_ambient_counts_before_after",
        caption="Per-cell total counts before and after ambient correction.",
        render=render,
    )
    return (
        "## Ambient RNA correction\n"
        "\n"
        + rho_note
        + counts_note
        + fig
        + "\n"
        f"{_md_table(['parameter', 'value'], rows)}\n"
    )


def _rna_section(
    run_dir: Path,
    params: dict[str, Any],
    counts: dict[str, Any],
    *,
    render: _FigRender = _DEFAULT_FIG_RENDER,
) -> str:
    from .run_paths import RunPaths
    s1 = RunPaths(run_dir).stage_dir("s1_rna_qc")
    pre = s1 / "qc_metrics_pre.parquet"
    post = s1 / "qc_metrics_post.parquet"
    if not (pre.exists() and post.exists()):
        return "## RNA QC\n\n_(artifacts not available)_\n"

    pre_df = pd.read_parquet(pre)
    post_df = pd.read_parquet(post)
    n_pre = len(pre_df)
    n_post = len(post_df)
    n_rm = n_pre - n_post

    thresholds = _md_table(
        ["parameter", "value"],
        [
            ["total_counts_min", _param(params, "s1_rna_qc.total_counts_min")],
            ["total_counts_max", _param(params, "s1_rna_qc.total_counts_max")],
            ["n_genes_min",      _param(params, "s1_rna_qc.n_genes_min")],
            ["n_genes_max",      _param(params, "s1_rna_qc.n_genes_max")],
            ["pct_counts_mt_max", _param(params, "s1_rna_qc.pct_counts_mt_max")],
            ["pct_counts_ribo_max", _param(params, "s1_rna_qc.pct_counts_ribo_max")],
            ["n_mt_genes_detected", _param(params, "s1_rna_qc.n_mt_genes_detected")],
            ["n_ribo_genes_detected", _param(params, "s1_rna_qc.n_ribo_genes_detected")],
        ],
    )

    stat_rows: list[list[Any]] = []
    for col in ("n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"):
        if col in post_df.columns:
            stat_rows.append(_stats_row(col, post_df[col].to_numpy()))
    stats = _md_table(["metric", "mean", "median", "min", "max"], stat_rows) if stat_rows else ""

    figs = (
        _fig_block(run_dir, "s1_rna_qc_violin_pre",
                   caption="RNA QC metrics before filtering.", render=render)
        + _fig_block(run_dir, "s1_rna_qc_violin_post",
                      caption="RNA QC metrics after filtering.", render=render)
    )
    return (
        "## RNA quality filtering\n"
        "\n"
        "Removes outliers and low-quality cells using MAD-based bounds on total UMI "
        "counts and detected genes, plus ceilings on mitochondrial (MT) and ribosomal "
        "read fractions.\n"
        "\n"
        f"- Cells before filtering: **{n_pre}**\n"
        f"- Cells retained:         **{n_post}**\n"
        f"- Removed:                **{n_rm}** ({_pct(n_rm, n_pre)})\n"
        + figs
        + "\n"
        "### Thresholds used\n\n"
        f"{thresholds}\n"
        "\n"
        "### Summary statistics (retained cells)\n\n"
        f"{stats}\n"
    )


def _atac_section(
    run_dir: Path,
    params: dict[str, Any],
    counts: dict[str, Any],
    *,
    render: _FigRender = _DEFAULT_FIG_RENDER,
) -> str:
    from .run_paths import RunPaths
    s2 = RunPaths(run_dir).stage_dir("s2_atac_qc")
    summary_json = s2 / "qc_summary.json"
    atac_h5ad = s2 / "atac_qc.h5ad"
    if not summary_json.exists():
        return "## ATAC QC\n\n_(artifacts not available)_\n"

    summary = json.loads(summary_json.read_text())
    n_pre = int(summary.get("n_cells_pre", 0))       # post-import, pre-S2-filter
    n_after_3m = summary.get("n_cells_after_3m_filter")  # after n_frag/TSS/NS, before FRiP
    n_post = int(summary.get("n_cells_post", 0))
    n_rm = n_pre - n_post
    peak_source = summary.get("peak_source")
    atac_raw = counts.get("atac_raw_barcodes")
    snap_drop = (atac_raw - n_pre) if (atac_raw is not None) else None

    frip_min_val = _param(params, "s2_atac_qc.frip_min")
    frip_threshold_note = frip_min_val if peak_source else f"{frip_min_val} _(not applied — no peaks available)_"
    thresholds = _md_table(
        ["parameter", "value"],
        [
            ["n_fragments_min",    _param(params, "s2_atac_qc.n_fragments_min")],
            ["n_fragments_max",    _param(params, "s2_atac_qc.n_fragments_max")],
            ["tss_enrichment_min", _param(params, "s2_atac_qc.tss_enrichment_min")],
            ["tss_enrichment_max", _param(params, "s2_atac_qc.tss_enrichment_max")],
            ["nucleosome_signal_max", _param(params, "s2_atac_qc.nucleosome_signal_max")],
            ["frip_min",           frip_threshold_note],
        ],
    )

    # Summary stats — read the post-QC SnapATAC2 AnnData obs
    stat_rows: list[list[Any]] = []
    warnings: list[str] = []
    if atac_h5ad.exists():
        try:
            import snapatac2 as snap
            adata = snap.read(str(atac_h5ad))
            obs = adata.obs[:].to_pandas()
            try: adata.close()
            except Exception: pass
            for src_col, label in [
                ("n_fragment", "fragment_count"),
                ("tsse", "tss_enrichment"),
                ("nucleosome_signal", "nucleosome_signal"),
                ("frip", "frip"),
            ]:
                if src_col in obs.columns:
                    stat_rows.append(_stats_row(label, obs[src_col].to_numpy()))
        except BaseException as e:
            warnings.append(f"_Could not read ATAC AnnData for summary stats: {e}_")
    stats = _md_table(["metric", "mean", "median", "min", "max"], stat_rows) if stat_rows else "_(no stats)_"
    warn_block = ("\n" + "\n".join(warnings) + "\n") if warnings else ""

    import_note = ""
    if snap_drop is not None and snap_drop > 0:
        import_note = (
            f"- **Fragment-count pre-filter (at import):** {snap_drop} barcodes removed "
            f"({atac_raw} → {n_pre}) — cells with too few fragments are dropped before "
            f"quality metrics are computed.\n"
        )

    # Two-stage waterfall when FRiP was computed
    if n_after_3m is not None:
        n_after_3m = int(n_after_3m)
        n_rm_3m = n_pre - n_after_3m
        n_rm_frip = n_after_3m - n_post
        count_lines = (
            f"- Cells before filtering:       **{n_pre}**\n"
            f"- After n_frag/TSS/NS filter:   **{n_after_3m}** (removed {n_rm_3m}, {_pct(n_rm_3m, n_pre)})\n"
            f"- After FRiP filter:            **{n_post}** (removed {n_rm_frip}, {_pct(n_rm_frip, n_after_3m)})\n"
        )
    else:
        count_lines = (
            f"- Cells before filtering: **{n_pre}**\n"
            f"- Cells retained:         **{n_post}**\n"
            f"- Removed:                **{n_rm}** ({_pct(n_rm, n_pre)})\n"
        )

    peak_note = f"\n_Peak source for FRiP: {peak_source}._\n" if peak_source else ""

    from .figures import FRIP_DISTRIBUTION_TITLE, TSS_PROFILE_CAPTION, TSS_PROFILE_TITLE
    fig_fsd_frip = _fig_pair_block(
        run_dir,
        ("s2_atac_qc_fragment_size_distribution",
         "ATAC fragment size distribution after QC filtering."),
        ("s2_atac_qc_frip_histogram",
         f"{FRIP_DISTRIBUTION_TITLE}; dashed line marks the filter threshold."),
        render=render,
    )
    fig_tss = _fig_block(
        run_dir, "s2_atac_qc_tss_enrichment_profile",
        caption=f"{TSS_PROFILE_TITLE}; {TSS_PROFILE_CAPTION}",
        render=render,
    )
    return (
        "## ATAC quality filtering\n"
        "\n"
        "Removes low-quality cells using MAD-based bounds on fragment counts, plus "
        "TSS enrichment, nucleosome-signal, and FRiP thresholds.\n"
        "\n"
        f"{import_note}"
        + count_lines
        + "\n"
        "### Thresholds used\n\n"
        f"{thresholds}\n"
        f"{peak_note}\n"
        "### Summary statistics (retained cells)\n\n"
        f"{stats}\n"
        f"{warn_block}"
        + fig_fsd_frip
        + fig_tss
    )


def _workflow_branch(run_dir: Path, params: dict[str, Any]) -> str:
    from . import provenance as _prov
    from .run_paths import RunPaths
    branch = _prov.get_value(str(RunPaths(run_dir).parameters_yaml), "plan.workflow_branch")
    if branch:
        return str(branch)
    return str(params.get("plan.workflow_branch", "paired"))


def _plan_dataset_assay_line(run_dir: Path) -> str:
    """One-line dataset type and assay from P1 context (same fields as plan review summary)."""
    from .run_paths import RunPaths
    ctx_path = RunPaths(run_dir).artifact("p1_context", "context_extraction.json")
    if not ctx_path.exists():
        return ""
    try:
        ctx = json.loads(ctx_path.read_text(encoding="utf-8"))
    except Exception:
        return ""
    fields = ctx.get("fields", {})
    dtype = (fields.get("modality_type") or {}).get("value") or "unknown"
    assay = (fields.get("assay_type") or {}).get("value") or "unknown"
    return f"**Dataset type:** {dtype} · **Assay:** {assay}"


def _qc_review_intro(run_name: str, run_dir: Path) -> str:
    lines = [
        "# QC review checkpoint",
        "",
        f"Review QC plots in `{run_name}/deliverables/checkpoint/qc_review`. "
        f"For a rendered report with images, open **qc_summary_{run_name}.html** in this folder "
        "(generated alongside this file). Approve when satisfied, or revise thresholds "
        "and re-run affected stages before approving.",
        "",
    ]
    context = _plan_dataset_assay_line(run_dir)
    if context:
        lines += [context, ""]
    return "\n".join(lines)


def _qc_review_actions(branch: str) -> str:
    lines = [
        "## How to approve or revise",
        "",
        "### Review QC filtering thresholds",
        "",
        "Look at the RNA and ATAC data distributions. Decide whether the current QC "
        "filtering thresholds look appropriate for the data, or whether they should be "
        "made stricter or more permissive.",
        "",
        "If you want to adjust the thresholds, tell the agent which stage to revise and how:",
        "",
        "- **S1 (RNA QC)** — UMI count, gene count, mitochondrial fraction, and ribosomal fraction bounds",
        "- **S2 (ATAC QC)** — fragment count, TSS enrichment, nucleosome signal, and FRiP",
        "- **S3 (doublets)** — RNA Scrublet and ATAC SnapATAC2 score cutoffs",
        "",
        "The agent can change the settings for the affected stage, re-run downstream steps "
        "as needed, and regenerate the reports.",
        "",
        "### Approve and continue to downstream analysis",
        "",
        "If the QC filters look acceptable, tell the agent to approve this checkpoint and "
        "continue to downstream dimensionality reduction and clustering.",
        "",
    ]
    return "\n".join(lines)


def _doublet_section(
    run_dir: Path,
    params: dict[str, Any],
    counts: dict[str, Any],
    *,
    render: _FigRender = _DEFAULT_FIG_RENDER,
) -> str:
    """Build a corrected doublet summary from the raw per-cell calls.parquet.

    Issues with the existing `overlap_summary.json`:
      - Applies fillna(False) before classifying, so "not evaluated" collapses
        into "neither flagged".
      - `n_removed` is the distinct-barcode count in the merged set, not the
        per-modality removal count.

    This section re-derives the overlap over cells evaluated by BOTH detectors
    (non-null scores in both columns) and reports per-modality removals by
    comparing against the post-doublet h5ads.
    """
    from .run_paths import RunPaths
    s3 = RunPaths(run_dir).stage_dir("s3_doublets")
    calls_path = s3 / "calls.parquet"
    overlap_path = s3 / "overlap_summary.json"
    if not (calls_path.exists() and overlap_path.exists()):
        return "## Doublet removal\n\n_(artifacts not available)_\n"

    calls = pd.read_parquet(calls_path)
    # Evaluation status per cell (based on whether each detector scored it)
    rna_scored = calls["scrublet_score"].notna()
    atac_scored = calls["atac_doublet_score"].notna()
    both_scored = rna_scored & atac_scored
    only_rna_scored = rna_scored & ~atac_scored
    only_atac_scored = ~rna_scored & atac_scored

    # Per-detector totals (across all cells the detector saw)
    n_rna_flag_total = int(calls.loc[rna_scored, "scrublet_is_doublet"].fillna(False).sum())
    n_atac_flag_total = int(calls.loc[atac_scored, "atac_is_doublet"].fillna(False).sum())

    # Corrected four-way overlap — restricted to cells scored by BOTH detectors
    cells_both = calls[both_scored]
    rna_flag_b = cells_both["scrublet_is_doublet"].fillna(False)
    atac_flag_b = cells_both["atac_is_doublet"].fillna(False)
    both_n = int(((rna_flag_b) & (atac_flag_b)).sum())
    rna_only_n = int(((rna_flag_b) & (~atac_flag_b)).sum())
    atac_only_n = int(((~rna_flag_b) & (atac_flag_b)).sum())
    neither_n = int(((~rna_flag_b) & (~atac_flag_b)).sum())
    n_both = int(both_scored.sum())

    # Per-modality doublet removals: compare S2/S1 post-QC cells to post-doublet cells
    n_removed_rna = None
    if counts.get("rna_qc_post") is not None and counts.get("rna_post_doublet") is not None:
        n_removed_rna = counts["rna_qc_post"] - counts["rna_post_doublet"]
    n_removed_atac = None
    if counts.get("atac_qc_post") is not None and counts.get("atac_post_doublet") is not None:
        n_removed_atac = counts["atac_qc_post"] - counts["atac_post_doublet"]
    # Distinct flagged barcodes (in the union merged set) — this is what the
    # raw overlap_summary.json previously called "n_removed".
    n_distinct_flagged = int(((calls["scrublet_is_doublet"].fillna(False)) |
                              (calls["atac_is_doublet"].fillna(False))).sum())

    overlap_summary = json.loads(overlap_path.read_text())
    policy = (_param(params, "s3_doublets.removal_policy")
              or overlap_summary.get("recommended_policy")
              or overlap_summary.get("chosen_policy")
              or overlap_summary.get("policy") or "unspecified")
    branch = _workflow_branch(run_dir, params)

    policy_note = ""

    if branch == "paired":
        policy_note = (
            "\n"
            "### Cross-modal policy (paired)\n"
            "\n"
            "- Applied policy: **union** — remove if either RNA (Scrublet) or ATAC "
            "(SnapATAC2) flags a doublet.\n"
            + f"- Joint cells after doublet removal: **{counts.get('n_cells_joint', 'n/a')}**\n"
        )
    elif branch == "separate":
        policy_note = (
            "\n"
            "### Per-modality removal (separate branch)\n"
            "\n"
            "RNA and ATAC doublet calls are applied independently; each modality "
            "keeps its own survivor set.\n"
        )
    else:
        policy_note = (
            "\n"
            f"### Single-modality removal ({branch} branch)\n"
            "\n"
            "Doublet filtering runs on the present modality only.\n"
        )

    rna_before = int(rna_scored.sum())
    atac_before = int(atac_scored.sum())

    rna_thr = (_param(params, "s3_doublets.rna_doublet_score_threshold") if params else None)
    if rna_thr is None and params is not None:
        rna_thr = 0.25
    atac_thr = None
    if params:
        atac_thr = (_param(params, "s3_doublets.atac_doublet_probability_threshold")
                    or _param(params, "s3_doublets.atac_doublet_threshold")
                    or _param(params, "s3_doublets.atac_doublet_score_threshold"))
        if atac_thr is None:
            atac_thr = 0.5
    threshold_lines = ""
    if rna_thr is not None or atac_thr is not None:
        thr_rows: list[list[Any]] = []
        if rna_thr is not None:
            thr_rows.append(["rna_doublet_score_threshold", rna_thr])
        if atac_thr is not None:
            thr_rows.append(["atac_doublet_probability_threshold", atac_thr])
        threshold_lines = (
            "\n"
            "### Doublet thresholds\n"
            "\n"
            f"{_md_table(['parameter', 'value'], thr_rows)}\n"
        )

    summary_table = ""
    if branch == "paired":
        both_scored_label = f"n={n_both}*"
        doublet_tbl = _md_table(
            ["", "RNA before", "ATAC before", "RNA-only flagged",
             "ATAC-only flagged", "Both flagged", "**Retained after removal**"],
            [[both_scored_label, rna_before, atac_before, rna_only_n, atac_only_n,
              both_n, neither_n]],
        )
        footnote = (
            "- *Cell barcodes that were evaluated for doublets by both the RNA (Scrublet) "
            "and ATAC (SnapATAC2) detectors. **Retained** = neither flagged by any detector.\n"
        )
        summary_table = (
            "\n"
            f"{doublet_tbl}\n"
            "\n"
            f"{footnote}"
        )
    else:
        summary_table = (
            "\n"
            f"- Applied removal policy: **{policy}**\n"
            f"- RNA evaluated: **{rna_before}** (flagged **{n_rna_flag_total}**)\n"
            f"- ATAC evaluated: **{atac_before}** (flagged **{n_atac_flag_total}**)\n"
        )

    doublet_figs = _fig_pair_block(
        run_dir,
        ("post_qc_review_doublet_rna", "RNA doublet scores (Scrublet)."),
        ("post_qc_review_doublet_atac", "ATAC doublet scores (SnapATAC2)."),
        render=render,
    )

    return (
        "## Doublet removal\n"
        f"{threshold_lines}"
        f"{policy_note}"
        f"{summary_table}"
        f"{doublet_figs}"
    )


def _final_section(run_dir: Path, counts: dict[str, Any]) -> str:
    rna_pd = counts.get("rna_post_doublet")
    atac_pd = counts.get("atac_post_doublet")
    rna_final = counts.get("rna_final")
    atac_final = counts.get("atac_final")

    if rna_final is None and atac_final is None:
        return "## Final retained dataset\n\n_(processed objects not available)_\n"

    from .run_paths import RunPaths
    s8 = RunPaths(run_dir).stage_dir("s8_umap")
    branch = "paired (processed.h5mu)" if (s8 / "processed.h5mu").exists() else "separate (two h5ads)"

    # Cluster counts
    try:
        import mudata as mu
        import anndata as ad
        if (s8 / "processed.h5mu").exists():
            m = mu.read_h5mu(str(s8 / "processed.h5mu"))
            rna = m.mod.get("rna")
            atac = m.mod.get("atac")
        else:
            rna = ad.read_h5ad(s8 / "rna_processed.h5ad")
            atac = ad.read_h5ad(s8 / "atac_processed.h5ad")
    except Exception as e:
        return f"## Final retained dataset\n\n_Error loading processed object: {e}_\n"

    def _n_clusters(ad_, label_col: str) -> int | None:
        if ad_ is None or label_col not in ad_.obs.columns:
            return None
        return int(ad_.obs[label_col].astype(str).nunique())

    rna_k = _n_clusters(rna, "leiden_rna")
    atac_k = _n_clusters(atac, "leiden_atac")

    matched = (
        rna is not None and atac is not None
        and rna.n_obs == atac.n_obs and rna.n_obs > 0
        and set(rna.obs_names) == set(atac.obs_names)
    )
    match_str = "yes (barcodes aligned)" if matched else "no"

    # Drop from S3 (post-doublet, post-intersection) to S8 (final).
    rna_lost = (rna_pd - rna_final) if (rna_pd is not None and rna_final is not None) else None
    atac_lost = (atac_pd - atac_final) if (atac_pd is not None and atac_final is not None) else None
    joint = counts.get("n_cells_joint")

    return (
        "## Final retained dataset (S8)\n"
        "\n"
        f"- Output: **{branch}**\n"
        "\n"
        "### Joint cell set (S3 paired intersection)\n"
        "\n"
        f"- Joint barcodes after S3 intersection: **{joint if joint is not None else 'n/a (non-paired branch)'}**\n"
        f"- RNA entering downstream stages:  **{rna_pd if rna_pd is not None else 'n/a'}**\n"
        f"- ATAC entering downstream stages: **{atac_pd if atac_pd is not None else 'n/a'}**\n"
        f"- Cells dropped between S3 and final assembly: RNA **{rna_lost if rna_lost is not None else 'n/a'}**, "
        f"ATAC **{atac_lost if atac_lost is not None else 'n/a'}** "
        "(should be zero on paired branch — S8 assembly is a safety no-op).\n"
        "\n"
        "### Final counts\n"
        "\n"
        f"- RNA cells:  **{rna_final if rna_final is not None else 'n/a'}**\n"
        f"- ATAC cells: **{atac_final if atac_final is not None else 'n/a'}**\n"
        f"- Modalities matched: **{match_str}**\n"
        f"- RNA clusters (leiden_rna):   **{rna_k if rna_k is not None else 'n/a'}** _(diagnostic, per-modality)_\n"
        f"- ATAC clusters (leiden_atac): **{atac_k if atac_k is not None else 'n/a'}** _(diagnostic, per-modality)_\n"
    )


# ---------------------------------------------------------------------------
# Top-level
# ---------------------------------------------------------------------------

def _qc_review_result_sections(
    run_dir: Path,
    params: dict[str, Any],
    counts: dict[str, Any],
    render: _FigRender,
    workflow_branch: str,
) -> list[str]:
    """QC metrics and figures only (no workflow intro or approve/revise instructions)."""
    return [
        _flow_section(run_dir, counts, workflow_branch, include_final_stage=False)
        + _flow_figures(run_dir, render=render),
        _ambient_section(run_dir, params, counts, render=render),
        _rna_section(run_dir, params, counts, render=render),
        _atac_section(run_dir, params, counts, render=render),
        _doublet_section(run_dir, params, counts, render=render),
    ]


def build_qc_review(run_dir: Path | str) -> str:
    """Markdown for the QC review user checkpoint (after S3, before S4/S5)."""
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    rp = RunPaths(run_dir)
    params_path = rp.parameters_yaml
    params = yaml.safe_load(params_path.read_text()) if params_path.exists() else {}
    counts = _stage_counts(run_dir)
    branch = _workflow_branch(run_dir, params)
    render = _FigRender(md_parent=rp.deliv_qc_review, embed_figures=True)

    sections = [
        _qc_review_intro(run_dir.name, run_dir),
        *_qc_review_result_sections(run_dir, params, counts, render, branch),
        _qc_review_actions(branch),
    ]
    return "\n".join(sections).rstrip() + "\n"


def build_qc_review_report(run_dir: Path | str) -> str:
    """Results-only markdown for HTML report (title + metrics/figures)."""
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    rp = RunPaths(run_dir)
    params_path = rp.parameters_yaml
    params = yaml.safe_load(params_path.read_text()) if params_path.exists() else {}
    counts = _stage_counts(run_dir)
    render = _FigRender(md_parent=rp.deliv_qc_review, embed_figures=True)
    branch = _workflow_branch(run_dir, params)
    title = f"# QC review — {run_dir.name}\n"
    context = _plan_dataset_assay_line(run_dir)
    header = title + (f"\n{context}\n" if context else "")
    sections = [header, *_qc_review_result_sections(run_dir, params, counts, render, branch)]
    return "\n".join(sections).rstrip() + "\n"


def _qc_html_styles() -> str:
    return (
        "body { font-family: system-ui, sans-serif; max-width: 56rem; "
        "margin: 2rem auto; padding: 0 1rem; line-height: 1.5; color: #111; }\n"
        "h1 { font-size: 1.5rem; margin-bottom: 1.25rem; }\n"
        "h2 { font-size: 1.2rem; margin-top: 2rem; border-bottom: 1px solid #e5e5e5; "
        "padding-bottom: 0.25rem; }\n"
        "h3 { font-size: 1.05rem; margin-top: 1.25rem; }\n"
        "table { border-collapse: collapse; margin: 0.75rem 0; width: 100%; "
        "font-size: 0.92rem; }\n"
        "th, td { border: 1px solid #ccc; padding: 0.35rem 0.6rem; text-align: left; }\n"
        "th { background: #f8f8f8; }\n"
        "ul { margin: 0.5rem 0 1rem 1.25rem; }\n"
        "p { margin: 0.5rem 0; }\n"
        "li.qc-footnote { font-size: 0.92rem; color: #444; }\n"
        ".qc-section { margin-bottom: 1.5rem; }\n"
        ".qc-side-by-side { display: flex; flex-wrap: wrap; gap: 1rem; "
        "align-items: flex-start; margin: 1rem 0; }\n"
        ".qc-side-by-side .qc-panel-table { flex: 1 1 260px; min-width: 0; }\n"
        ".qc-side-by-side .qc-panel-plot { flex: 0 1 380px; max-width: 52%; "
        "min-width: 240px; }\n"
        ".qc-plot-pair { display: flex; flex-wrap: wrap; gap: 1rem; margin: 1rem 0; "
        "align-items: stretch; }\n"
        ".qc-plot-pair .qc-figure { flex: 1 1 45%; min-width: 260px; margin: 0; "
        "display: flex; flex-direction: column; }\n"
        ".qc-plot-pair .qc-figure img { width: 100%; height: 300px; object-fit: contain; "
        "object-position: center top; display: block; }\n"
        ".qc-plots-below { margin-top: 1rem; }\n"
        ".qc-figure { margin: 0; }\n"
        ".qc-figure img { width: 100%; height: auto; display: block; }\n"
        ".qc-figure figcaption { font-size: 0.85rem; color: #555; margin-top: 0.35rem; }\n"
    )


def _html_figure(fig_dir: Path, stem: str, caption: str) -> str:
    png = fig_dir / f"{stem}.png"
    if not png.exists():
        return ""
    data = base64.standard_b64encode(png.read_bytes()).decode("ascii")
    alt = html_module.escape(caption)
    return (
        f'<figure class="qc-figure"><img src="data:image/png;base64,{data}" '
        f'alt="{alt}"><figcaption>{alt}</figcaption></figure>'
    )


def _md_html(md: str) -> str:
    return _markdown_to_html(md.strip()) if md.strip() else ""


def _html_flow_section(
    run_dir: Path, counts: dict[str, Any], workflow_branch: str, fig_dir: Path,
) -> str:
    table_md = _flow_section(run_dir, counts, workflow_branch, include_final_stage=False)
    plot = _html_figure(fig_dir, "post_qc_review_cell_counts",
                        "Cell counts across preprocessing (RNA and ATAC).")
    plots = f'<div class="qc-plots-below">{plot}</div>' if plot else ""
    return (
        '<section class="qc-section qc-flow">'
        f"{_md_html(table_md)}{plots}</section>"
    )


def _html_ambient_section(
    run_dir: Path, params: dict[str, Any], counts: dict[str, Any], fig_dir: Path,
) -> str:
    md = _ambient_section(run_dir, params, counts, render=_FigRender(embed_figures=False))
    if "## Ambient RNA correction\n" not in md or "| parameter | value |" not in md:
        return f'<section class="qc-section qc-ambient">{_md_html(md)}</section>'

    parts = md.split("## Ambient RNA correction\n", 1)[1]
    table_marker = "| parameter | value |"
    intro, _, rest = parts.partition(table_marker)
    table_lines = "| parameter | value |" + rest.split("\n\n", 1)[0]
    intro = intro.strip()
    plot = _html_figure(
        fig_dir, "s1a_ambient_counts_before_after",
        "Per-cell total counts before and after ambient correction.",
    )
    return (
        '<section class="qc-section qc-ambient">'
        "<h2>Ambient RNA correction</h2>"
        f"{_md_html(intro)}"
        '<div class="qc-side-by-side">'
        f'<div class="qc-panel-table">{_md_html(table_lines)}</div>'
        f'<div class="qc-panel-plot">{plot}</div>'
        "</div></section>"
    )


def _html_rna_section(
    run_dir: Path, params: dict[str, Any], counts: dict[str, Any], fig_dir: Path,
) -> str:
    from .run_paths import RunPaths
    s1 = RunPaths(run_dir).stage_dir("s1_rna_qc")
    pre = s1 / "qc_metrics_pre.parquet"
    post = s1 / "qc_metrics_post.parquet"
    if not (pre.exists() and post.exists()):
        return _md_html(_rna_section(run_dir, params, counts, render=_FigRender(embed_figures=False)))

    pre_df = pd.read_parquet(pre)
    post_df = pd.read_parquet(post)
    n_pre, n_post = len(pre_df), len(post_df)
    n_rm = n_pre - n_post

    thresholds = _md_table(
        ["parameter", "value"],
        [
            ["total_counts_min", _param(params, "s1_rna_qc.total_counts_min")],
            ["total_counts_max", _param(params, "s1_rna_qc.total_counts_max")],
            ["n_genes_min", _param(params, "s1_rna_qc.n_genes_min")],
            ["n_genes_max", _param(params, "s1_rna_qc.n_genes_max")],
            ["pct_counts_mt_max", _param(params, "s1_rna_qc.pct_counts_mt_max")],
            ["pct_counts_ribo_max", _param(params, "s1_rna_qc.pct_counts_ribo_max")],
            ["n_mt_genes_detected", _param(params, "s1_rna_qc.n_mt_genes_detected")],
            ["n_ribo_genes_detected", _param(params, "s1_rna_qc.n_ribo_genes_detected")],
        ],
    )
    stat_rows: list[list[Any]] = []
    for col in ("n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"):
        if col in post_df.columns:
            stat_rows.append(_stats_row(col, post_df[col].to_numpy()))
    stats = _md_table(["metric", "mean", "median", "min", "max"], stat_rows) if stat_rows else ""

    plots = (
        _html_figure(fig_dir, "s1_rna_qc_violin_pre", "RNA QC metrics before filtering.")
        + _html_figure(fig_dir, "s1_rna_qc_violin_post", "RNA QC metrics after filtering.")
    )
    intro = (
        "Removes outliers and low-quality cells using MAD-based bounds on total UMI "
        "counts and detected genes, plus ceilings on mitochondrial (MT) and ribosomal "
        "read fractions.\n\n"
        f"- Cells before filtering: **{n_pre}**\n"
        f"- Cells retained:         **{n_post}**\n"
        f"- Removed:                **{n_rm}** ({_pct(n_rm, n_pre)})\n"
    )
    return (
        '<section class="qc-section qc-rna">'
        "<h2>RNA quality filtering</h2>"
        f"{_md_html(intro)}"
        "<h3>Thresholds used</h3>"
        f"{_md_html(thresholds)}"
        "<h3>Summary statistics (retained cells)</h3>"
        f"{_md_html(stats)}"
        f'<div class="qc-plots-below">{plots}</div>'
        "</section>"
    )


def _html_atac_section(
    run_dir: Path, params: dict[str, Any], counts: dict[str, Any], fig_dir: Path,
) -> str:
    from .run_paths import RunPaths
    s2 = RunPaths(run_dir).stage_dir("s2_atac_qc")
    summary_json = s2 / "qc_summary.json"
    if not summary_json.exists():
        return _md_html(_atac_section(run_dir, params, counts, render=_FigRender(embed_figures=False)))

    summary = json.loads(summary_json.read_text())
    n_pre = int(summary.get("n_cells_pre", 0))
    n_post = int(summary.get("n_cells_post", 0))
    n_rm = n_pre - n_post
    atac_raw = counts.get("atac_raw_barcodes")
    snap_drop = (atac_raw - n_pre) if (atac_raw is not None) else None

    peak_source = summary.get("peak_source")
    frip_min_val = _param(params, "s2_atac_qc.frip_min")
    frip_threshold_note = frip_min_val if peak_source else f"{frip_min_val} _(not applied — no peaks available)_"
    thresholds = _md_table(
        ["parameter", "value"],
        [
            ["n_fragments_min", _param(params, "s2_atac_qc.n_fragments_min")],
            ["n_fragments_max", _param(params, "s2_atac_qc.n_fragments_max")],
            ["tss_enrichment_min", _param(params, "s2_atac_qc.tss_enrichment_min")],
            ["tss_enrichment_max", _param(params, "s2_atac_qc.tss_enrichment_max")],
            ["nucleosome_signal_max", _param(params, "s2_atac_qc.nucleosome_signal_max")],
            ["frip_min", frip_threshold_note],
        ],
    )
    stat_rows: list[list[Any]] = []
    atac_h5ad = s2 / "atac_qc.h5ad"
    if atac_h5ad.exists():
        try:
            import snapatac2 as snap
            adata = snap.read(str(atac_h5ad))
            obs = adata.obs[:].to_pandas()
            try:
                adata.close()
            except Exception:
                pass
            for src_col, label in [
                ("n_fragment", "fragment_count"),
                ("tsse", "tss_enrichment"),
                ("nucleosome_signal", "nucleosome_signal"),
                ("frip", "frip"),
            ]:
                if src_col in obs.columns:
                    stat_rows.append(_stats_row(label, obs[src_col].to_numpy()))
        except BaseException:
            pass
    stats = _md_table(["metric", "mean", "median", "min", "max"], stat_rows) if stat_rows else "_(no stats)_"

    import_note = ""
    if snap_drop is not None and snap_drop > 0:
        import_note = (
            f"- **Fragment-count pre-filter (at import):** {snap_drop} barcodes removed "
            f"({atac_raw} → {n_pre}) — cells with too few fragments are dropped before "
            f"quality metrics are computed.\n"
        )
    intro = (
        "Removes low-quality cells using MAD-based bounds on fragment counts, plus "
        "TSS enrichment, nucleosome-signal, and FRiP thresholds.\n\n"
        f"{import_note}"
        f"- Cells before filtering: **{n_pre}**\n"
        f"- Cells retained:         **{n_post}**\n"
        f"- Removed:                **{n_rm}** ({_pct(n_rm, n_pre)})\n"
    )
    from .figures import FRIP_DISTRIBUTION_TITLE, TSS_PROFILE_CAPTION, TSS_PROFILE_TITLE
    plot_fsd = _html_figure(
        fig_dir, "s2_atac_qc_fragment_size_distribution",
        "ATAC fragment size distribution after QC filtering.",
    )
    plot_tss = _html_figure(
        fig_dir, "s2_atac_qc_tss_enrichment_profile",
        f"{TSS_PROFILE_TITLE}; {TSS_PROFILE_CAPTION}",
    )
    plot_frip = _html_figure(
        fig_dir, "s2_atac_qc_frip_histogram",
        f"{FRIP_DISTRIBUTION_TITLE}; dashed line marks the filter threshold.",
    )
    fsd_frip_pair = "".join(p for p in (plot_fsd, plot_frip) if p)
    pair_row = f'<div class="qc-plot-pair">{fsd_frip_pair}</div>' if fsd_frip_pair else ""
    tss_row = f'<div class="qc-plots-below">{plot_tss}</div>' if plot_tss else ""
    peak_note = f"\n_Peak source for FRiP: {peak_source}._\n" if peak_source else ""
    return (
        '<section class="qc-section qc-atac">'
        "<h2>ATAC quality filtering</h2>"
        f"{_md_html(intro)}"
        "<h3>Thresholds used</h3>"
        f"{_md_html(thresholds + peak_note)}"
        "<h3>Summary statistics (retained cells)</h3>"
        f"{_md_html(stats)}"
        f"{pair_row}{tss_row}"
        "</section>"
    )


def _html_doublet_section(
    run_dir: Path,
    params: dict[str, Any],
    counts: dict[str, Any],
    fig_dir: Path,
    workflow_branch: str,
) -> str:
    md = _doublet_section(run_dir, params, counts, render=_FigRender(embed_figures=False))
    policy_html = _md_html(md)

    rna_plot = _html_figure(
        fig_dir, "post_qc_review_doublet_rna", "RNA doublet scores (Scrublet).",
    )
    atac_plot = _html_figure(
        fig_dir, "post_qc_review_doublet_atac", "ATAC doublet scores (SnapATAC2).",
    )
    if workflow_branch == "paired":
        plots = "".join(p for p in (rna_plot, atac_plot) if p)
        plot_row = f'<div class="qc-plot-pair">{plots}</div>' if plots else ""
    else:
        stacked = "".join(p for p in (rna_plot, atac_plot) if p)
        plot_row = f'<div class="qc-plots-below">{stacked}</div>' if stacked else ""

    return f'<section class="qc-section qc-doublets">{policy_html}{plot_row}</section>'


def build_qc_review_html_body(run_dir: Path | str) -> str:
    """HTML report with layout: tables before plots; selected side-by-side panels."""
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    rp = RunPaths(run_dir)
    fig_dir = rp.deliv_qc_review
    params_path = rp.parameters_yaml
    params = yaml.safe_load(params_path.read_text()) if params_path.exists() else {}
    counts = _stage_counts(run_dir)
    branch = _workflow_branch(run_dir, params)

    title = f"<h1>QC summary — {html_module.escape(run_dir.name)}</h1>"
    context = _plan_dataset_assay_line(run_dir)
    context_html = f"<p>{_inline_md(context)}</p>" if context else ""

    parts = [
        title,
        context_html,
        _html_flow_section(run_dir, counts, branch, fig_dir),
        _html_ambient_section(run_dir, params, counts, fig_dir),
        _html_rna_section(run_dir, params, counts, fig_dir),
        _html_atac_section(run_dir, params, counts, fig_dir),
        _html_doublet_section(run_dir, params, counts, fig_dir, branch),
    ]
    return "\n".join(p for p in parts if p)


def write_qc_review_html(run_dir: Path | str) -> Path | None:
    """Write qc_summary_<run>.html — formatted results with embedded figures and layout CSS."""
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    rp = RunPaths(run_dir)
    out = rp.qc_summary_html
    body = build_qc_review_html_body(run_dir)
    doc = (
        "<!DOCTYPE html>\n<html lang=\"en\"><head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        f"<title>QC summary — {html_module.escape(run_dir.name)}</title>\n"
        "<style>\n"
        f"{_qc_html_styles()}"
        "</style>\n</head><body>\n"
        f"{body}\n</body></html>\n"
    )
    out.write_text(doc, encoding="utf-8")
    return out


def write_qc_review(run_dir: Path | str) -> Path:
    """Write editable qc_review_<run>.md and browser-friendly qc_summary_<run>.html."""
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    rp = RunPaths(run_dir)
    rp.deliv_qc_review.mkdir(parents=True, exist_ok=True)
    for legacy in (rp.deliv_qc_review / "qc_review.md", rp.deliv_qc_review / "qc_summary.html"):
        if legacy.exists():
            legacy.unlink()
    rp.qc_review_summary_md.write_text(build_qc_review(run_dir), encoding="utf-8")
    write_qc_review_html(run_dir)
    return rp.qc_review_summary_md


def build(run_dir: Path | str) -> str:
    """Full end-to-end QC summary (written at manifest to post_run/)."""
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    rp = RunPaths(run_dir)
    params_path = rp.parameters_yaml
    params = yaml.safe_load(params_path.read_text()) if params_path.exists() else {}
    counts = _stage_counts(run_dir)
    render = _FigRender(md_parent=rp.deliv_post_run)

    branch = _workflow_branch(run_dir, params)
    sections = [
        "# QC Summary",
        "",
        _flow_section(run_dir, counts, branch) + _flow_figures(run_dir, render=render),
        _ambient_section(run_dir, params, counts, render=render),
        _rna_section(run_dir, params, counts, render=render),
        _atac_section(run_dir, params, counts, render=render),
        _doublet_section(run_dir, params, counts, render=render),
        _final_section(run_dir, counts),
    ]
    return "\n".join(sections).rstrip() + "\n"


def write(run_dir: Path | str) -> Path:
    """Write the QC summary markdown directly to its canonical user-facing location."""
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    out = RunPaths(run_dir).qc_summary_md
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build(run_dir))
    return out
