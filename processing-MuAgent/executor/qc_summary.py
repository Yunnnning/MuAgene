"""Concise user-facing QC summary (markdown).

Audit-driven rewrite — reporting-only changes, no pipeline-logic edits:
    1. Cell-count flow table across stages (makes all transitions visible).
    2. Hidden SnapATAC2 import-stage cell drop is surfaced explicitly.
    3. Doublet overlap table restricted to cells evaluated by BOTH detectors;
       "not evaluated" is reported separately and no longer conflated with
       "not flagged".
    4. Per-modality doublet removal counts (computed from post-doublet h5ads).
    5. S8 paired-intersection loss shown in the "Final retained" section.
    6. Baselines relabelled from ambiguous "Cells before filtering" to
       stage-aware phrasing (e.g. "Cells entering this stage").
    7. Thresholds rounded to 2 decimal places; integer-typed values are
       displayed as integers.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


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


def _md_table(header: list[str], rows: list[list[Any]]) -> str:
    align = "|" + "|".join("---" for _ in header) + "|"
    h = "| " + " | ".join(header) + " |"
    body = "\n".join("| " + " | ".join(_fmt(x) for x in r) + " |" for r in rows)
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


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------

def _flow_section(counts: dict[str, Any]) -> str:
    """Cell-count flow across stages; makes transitions visible end-to-end."""
    def fmt(v):
        return "n/a" if v is None else str(int(v))

    rna_raw = counts["rna_raw"]
    atac_raw = counts["atac_raw_barcodes"]
    atac_after_import = counts["atac_after_snap_import"]
    # Derived gaps
    snap_drop = (atac_raw - atac_after_import) if (atac_raw is not None and atac_after_import is not None) else None

    rna_inter = (
        (counts["rna_post_doublet"] - counts["rna_final"])
        if (counts["rna_post_doublet"] is not None and counts["rna_final"] is not None) else None
    )
    atac_inter = (
        (counts["atac_post_doublet"] - counts["atac_final"])
        if (counts["atac_post_doublet"] is not None and counts["atac_final"] is not None) else None
    )

    rows = [
        ["1. raw (Cell Ranger / fragments)", fmt(rna_raw), fmt(atac_raw), "— / —"],
        ["2. after S0 ingest (paired intersection)",
         fmt(counts["rna_ingest"]), fmt(atac_raw),
         "RNA: filtered matrix cells intersected with ATAC barcodes"],
        ["3. after S1a ambient correction",
         fmt(counts["rna_after_ambient"]), fmt(atac_raw),
         "RNA: DecontX/SoupX correction (cells preserved; per-cell counts adjusted)"],
        ["4. after SnapATAC2 `import_fragments`",
         fmt(counts["rna_after_ambient"] if counts["rna_after_ambient"] is not None else counts["rna_ingest"]),
         fmt(atac_after_import),
         f"ATAC: SnapATAC2 `min_num_fragments` default drops {fmt(snap_drop)} cells"
         if snap_drop else "ATAC: SnapATAC2 import"],
        ["5. after S1 RNA QC  /  S2 ATAC QC",
         fmt(counts["rna_qc_post"]), fmt(counts["atac_qc_post"]),
         "explicit QC filters (RNA: MAD on counts/genes/mt + ribo ceiling; ATAC: TSS + n_frag + nucleosome_signal)"],
        ["6. after S3 doublet removal",
         fmt(counts["rna_post_doublet"]), fmt(counts["atac_post_doublet"]),
         "union / intersection policy, applied per modality"],
        ["7. after S8 paired-intersection (final)",
         fmt(counts["rna_final"]), fmt(counts["atac_final"]),
         f"RNA lost {fmt(rna_inter)}, ATAC lost {fmt(atac_inter)} at pairing"
         if rna_inter is not None else "paired MuData intersection"],
    ]
    return (
        "## Cell-count flow across stages\n\n"
        "Each row's RNA and ATAC columns describe the cell count _entering the "
        "next stage_ (equivalently, the count _after_ the stage named on that row).\n\n"
        f"{_md_table(['stage', 'RNA', 'ATAC', 'note'], rows)}\n"
    )


def _ambient_section(run_dir: Path, params: dict[str, Any], counts: dict[str, Any]) -> str:
    from .run_paths import RunPaths
    s1a = RunPaths(run_dir).stage_dir("s1a_ambient")
    summary_p = s1a / "summary.json"
    contam_p = s1a / "contamination.parquet"

    method = _param(params, "s1a_ambient.method")
    if method in (None, "none", "skipped_empty", "skipped_no_r"):
        note = {
            None: "_(stage did not run; legacy run or RNA absent)_",
            "none": "_user-disabled (method=none); no ambient correction applied._",
            "skipped_empty": "_RNA AnnData empty (atac_only branch); pass-through._",
            "skipped_no_r": ("_R / DecontX / SoupX not available at runtime; "
                              "stage degraded to pass-through. Install R + "
                              "BiocManager::install(c('celda','SoupX')) to enable._"),
        }[method]
        return "## Ambient RNA correction (S1a)\n\n" + note + "\n"

    summary: dict[str, Any] = {}
    if summary_p.exists():
        try:
            summary = json.loads(summary_p.read_text())
        except Exception:
            summary = {}

    rows = [
        ["method", method],
        ["median contamination", _param(params, "s1a_ambient.median_contamination")],
        ["high-contamination cells (rho>0.20)",
         _param(params, "s1a_ambient.n_high_contamination_cells")],
        ["max-contamination cap",
         summary.get("max_contam_cap", _param(params, "s1a_ambient.max_contamination"))],
        ["pre-correction total counts (sum)", summary.get("total_counts_pre")],
        ["post-correction total counts (sum)", summary.get("total_counts_post")],
    ]

    contam_stats = ""
    if contam_p.exists():
        try:
            df = pd.read_parquet(contam_p)
            v = df["contamination"].to_numpy()
            contam_stats = _md_table(
                ["metric", "mean", "median", "min", "max"],
                [_stats_row("contamination", v)],
            )
        except Exception:
            pass

    return (
        "## Ambient RNA correction (S1a)\n"
        "\n"
        "Decontaminated counts overwrite `.X` and `.layers['counts']`; the "
        "original counts are preserved in `.layers['counts_raw']`. Per-cell "
        "rho is in `.obs['ambient_contamination']` and `contamination.parquet`.\n"
        "\n"
        f"{_md_table(['parameter', 'value'], rows)}\n"
        "\n"
        + (("### Per-cell contamination distribution\n\n" + contam_stats + "\n") if contam_stats else "")
    )


def _rna_section(run_dir: Path, params: dict[str, Any], counts: dict[str, Any]) -> str:
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

    return (
        "## RNA QC (S1)\n"
        "\n"
        f"- Cells entering S1 (post-ingest, paired intersection): **{n_pre}**\n"
        f"- Cells retained after S1:                             **{n_post}**\n"
        f"- Removed by S1 filters:                               **{n_rm}** ({_pct(n_rm, n_pre)})\n"
        "\n"
        "### Thresholds used\n\n"
        f"{thresholds}\n"
        "\n"
        "### Summary statistics (post-S1, RNA-only cells)\n\n"
        f"{stats}\n"
    )


def _atac_section(run_dir: Path, params: dict[str, Any], counts: dict[str, Any]) -> str:
    from .run_paths import RunPaths
    s2 = RunPaths(run_dir).stage_dir("s2_atac_qc")
    summary_json = s2 / "qc_summary.json"
    atac_h5ad = s2 / "atac_qc.h5ad"
    if not summary_json.exists():
        return "## ATAC QC\n\n_(artifacts not available)_\n"

    summary = json.loads(summary_json.read_text())
    n_pre = int(summary.get("n_cells_pre", 0))       # post-import, pre-S2-filter
    n_post = int(summary.get("n_cells_post", 0))
    n_rm = n_pre - n_post
    atac_raw = counts.get("atac_raw_barcodes")
    snap_drop = (atac_raw - n_pre) if (atac_raw is not None) else None

    thresholds = _md_table(
        ["parameter", "value"],
        [
            ["n_fragments_min",    _param(params, "s2_atac_qc.n_fragments_min")],
            ["n_fragments_max",    _param(params, "s2_atac_qc.n_fragments_max")],
            ["tss_enrichment_min", _param(params, "s2_atac_qc.tss_enrichment_min")],
            ["nucleosome_signal_max", _param(params, "s2_atac_qc.nucleosome_signal_max")],
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
            ]:
                if src_col in obs.columns:
                    stat_rows.append(_stats_row(label, obs[src_col].to_numpy()))
        except Exception as e:
            warnings.append(f"_Could not read ATAC AnnData for summary stats: {e}_")
    stats = _md_table(["metric", "mean", "median", "min", "max"], stat_rows) if stat_rows else "_(no stats)_"
    warn_block = ("\n" + "\n".join(warnings) + "\n") if warnings else ""

    import_note = ""
    if snap_drop is not None and snap_drop > 0:
        import_note = (
            f"- **SnapATAC2 import-stage drop (hidden pre-filter):** "
            f"{snap_drop} cells ({atac_raw} → {n_pre}), due to SnapATAC2's "
            f"`snap.pp.import_fragments` default `min_num_fragments`.\n"
        )

    return (
        "## ATAC QC (S2)\n"
        "\n"
        f"{import_note}"
        f"- Cells entering S2 (post-SnapATAC2 import): **{n_pre}**\n"
        f"- Cells retained after S2:                   **{n_post}**\n"
        f"- Removed by S2 filters:                     **{n_rm}** ({_pct(n_rm, n_pre)})\n"
        "\n"
        "### Thresholds used\n\n"
        f"{thresholds}\n"
        "\n"
        "### Summary statistics (post-S2, ATAC cells)\n\n"
        f"{stats}\n"
        f"{warn_block}"
    )


def _doublet_section(run_dir: Path, counts: dict[str, Any]) -> str:
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
        return "## Doublets (S3)\n\n_(artifacts not available)_\n"

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
    policy = (overlap_summary.get("recommended_policy")
              or overlap_summary.get("chosen_policy") or "unspecified")

    return (
        "## Doublets (S3)\n"
        "\n"
        "### Per-detector flagged counts (evaluated cells only)\n"
        "\n"
        f"- Flagged by RNA (Scrublet), of {int(rna_scored.sum())} RNA-evaluated cells: **{n_rna_flag_total}**\n"
        f"- Flagged by ATAC (SnapATAC2), of {int(atac_scored.sum())} ATAC-evaluated cells: **{n_atac_flag_total}**\n"
        "\n"
        f"### Four-way overlap  (cells scored by BOTH detectors only — n={n_both})\n"
        "\n"
        f"{_md_table(['class', 'count'], [['RNA-only flagged', rna_only_n], ['ATAC-only flagged', atac_only_n], ['both flagged', both_n], ['neither flagged', neither_n]])}\n"
        "\n"
        "### Cells scored by only one detector (excluded from the overlap classification)\n"
        "\n"
        f"{_md_table(['status', 'count'], [['only RNA-scored (filtered at ATAC QC)', int(only_rna_scored.sum())], ['only ATAC-scored (filtered at RNA QC)', int(only_atac_scored.sum())]])}\n"
        "\n"
        "### Removal\n"
        "\n"
        f"- Policy: **{policy}**\n"
        f"- Removed from RNA (S3 RNA: {counts.get('rna_qc_post')} → {counts.get('rna_post_doublet')}): **{n_removed_rna if n_removed_rna is not None else 'n/a'}**\n"
        f"- Removed from ATAC (S3 ATAC: {counts.get('atac_qc_post')} → {counts.get('atac_post_doublet')}): **{n_removed_atac if n_removed_atac is not None else 'n/a'}**\n"
        f"- Distinct barcodes flagged across merged union set: **{n_distinct_flagged}** "
        f"_(for reference; this is the count the raw `overlap_summary.json` reported as `n_removed` and conflated with per-modality counts)_\n"
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

    # Paired-intersection delta
    rna_lost = (rna_pd - rna_final) if (rna_pd is not None and rna_final is not None) else None
    atac_lost = (atac_pd - atac_final) if (atac_pd is not None and atac_final is not None) else None

    return (
        "## Final retained dataset (S8)\n"
        "\n"
        f"- Output: **{branch}**\n"
        "\n"
        "### Paired-intersection step (S8)\n"
        "\n"
        f"- RNA entering S8 (post-doublet):  **{rna_pd if rna_pd is not None else 'n/a'}**\n"
        f"- ATAC entering S8 (post-doublet): **{atac_pd if atac_pd is not None else 'n/a'}**\n"
        f"- Cells dropped at pairing intersection: RNA **{rna_lost if rna_lost is not None else 'n/a'}**, "
        f"ATAC **{atac_lost if atac_lost is not None else 'n/a'}** "
        "(barcodes present in one modality's post-doublet object but not the other)\n"
        "\n"
        "### Final counts\n"
        "\n"
        f"- RNA cells:  **{rna_final if rna_final is not None else 'n/a'}**\n"
        f"- ATAC cells: **{atac_final if atac_final is not None else 'n/a'}**\n"
        f"- Modalities matched: **{match_str}**\n"
        f"- RNA clusters (leiden_rna):   **{rna_k if rna_k is not None else 'n/a'}**\n"
        f"- ATAC clusters (leiden_atac): **{atac_k if atac_k is not None else 'n/a'}**\n"
    )


# ---------------------------------------------------------------------------
# Top-level
# ---------------------------------------------------------------------------

def build(run_dir: Path | str) -> str:
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    params_path = RunPaths(run_dir).parameters_yaml
    params = yaml.safe_load(params_path.read_text()) if params_path.exists() else {}
    counts = _stage_counts(run_dir)

    sections = [
        "# QC Summary",
        "",
        _flow_section(counts),
        _ambient_section(run_dir, params, counts),
        _rna_section(run_dir, params, counts),
        _atac_section(run_dir, params, counts),
        _doublet_section(run_dir, counts),
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
