"""S1 — RNA QC thresholds (MAD-based) + filtering."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import scanpy as sc

from ..methods import mad_thresholds as _mad
from ..methods.qc_filter_stats import exclusive_removals
from .. import io as _io
from .. import provenance as _prov
from ..log import log_event


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s1_rna_qc"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"

    # S1a (ambient correction) is the canonical upstream for S1; fall back to
    # S0 directly only if S1a's artifact is missing for some reason (legacy
    # runs, or a hand-resumed pipeline). Both files contain raw integer counts
    # in `.layers["counts"]` and `.X`.
    s1a_path = run_dir / "internal" / "artifacts" / "s1a_ambient" / "rna_decontaminated.h5ad"
    s0_path = run_dir / "internal" / "artifacts" / "s0_ingest" / "rna_ingest.h5ad"
    branch = _prov.current_branch(str(params_path))
    if branch in ("paired", "separate", "rna_only") and not s1a_path.exists():
        raise FileNotFoundError(
            f"S1 RNA QC expected upstream S1a artifact at {s1a_path} but it is missing. "
            "Refusing to fall back to S0 ingest — that would skip ambient correction."
        )
    in_path = s1a_path if s1a_path.exists() else s0_path
    a = ad.read_h5ad(in_path)

    # Mito + ribosomal flagging.
    #
    # Symbols vary by species/case (Mt-/mt-/MT-/mt: for fly, etc.) and by
    # whether `var_names` are gene IDs (Ensembl) or symbols. Use a permissive
    # case-insensitive prefix set, then sanity-check by gene-count: a typical
    # mammalian mt gene set is ~13 genes; a typical ribosomal protein set is
    # ~80–90 genes. If we get suspiciously few mt genes we log a warning.
    mt_pat = r"(?i)^(mt[-:_]|mito[-:_]?)"      # mt-, MT-, Mt:, mito-, etc.
    ribo_pat = r"(?i)^(rps|rpl|mrps|mrpl)"     # cytoplasmic + mito ribosomal proteins
    a.var["mt"] = a.var_names.astype(str).str.contains(mt_pat, regex=True, na=False)
    a.var["ribo"] = a.var_names.astype(str).str.contains(ribo_pat, regex=True, na=False)
    n_mt = int(a.var["mt"].sum())
    n_ribo = int(a.var["ribo"].sum())
    if n_mt == 0:
        log_event(run_dir, {"stage": "s1_rna_qc", "event": "no_mt_genes_detected",
                             "note": "var_names may be Ensembl IDs; pct_counts_mt will be 0 and the MT filter is effectively disabled"})
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt", "ribo"], percent_top=None,
                                log1p=False, inplace=True)

    # Parameters from plan
    params = plan["stages"]["s1_rna_qc"]["parameters"]
    k_mad = params["k_mad"]["value"]
    pct_mt_k = params["pct_mt_k"]["value"]
    pct_mt_ceil = params["pct_mt_ceiling"]["value"]
    pct_mt_floor = params["pct_mt_floor"]["value"]
    min_cells = int(params["min_cells_per_gene"]["value"])
    min_counts_floor = params["min_counts_floor"]["value"]
    min_genes_floor = float(params.get("min_genes_floor", {}).get("value", 200))
    # Ribosomal upper-bound is recommended (not enforced strictly): some
    # tissues legitimately have very high ribo-protein expression.
    pct_ribo_max = float(params.get("pct_ribo_max", {}).get("value", 50.0))

    # Apply floor before MAD
    keep_floor = a.obs["total_counts"] >= min_counts_floor
    a_for_mad = a[keep_floor].copy()

    # Derive thresholds; absolute floors win when MAD lower bounds fall below them.
    c_lo, c_hi = _mad.log_mad_bounds(a_for_mad.obs["total_counts"].to_numpy(), k=k_mad)
    c_lo = max(c_lo, float(min_counts_floor))
    g_lo, g_hi = _mad.log_mad_bounds(a_for_mad.obs["n_genes_by_counts"].to_numpy(), k=k_mad)
    g_lo = max(g_lo, min_genes_floor)
    pct_mt_upper = _mad.upper_bound(a_for_mad.obs["pct_counts_mt"].to_numpy(),
                                     k=pct_mt_k, floor=pct_mt_floor, ceiling=pct_mt_ceil)

    # Record as provenance
    for key, value, rat in [
        ("s1_rna_qc.total_counts_min", float(c_lo),
         f"max(MAD lower bound on log1p(total_counts), min_counts_floor={min_counts_floor})"),
        ("s1_rna_qc.total_counts_max", float(c_hi), "MAD upper bound on log1p(total_counts)"),
        ("s1_rna_qc.n_genes_min", float(g_lo),
         f"max(MAD lower bound on log1p(n_genes), min_genes_floor={min_genes_floor})"),
        ("s1_rna_qc.n_genes_max", float(g_hi), "MAD upper bound on log1p(n_genes)"),
        ("s1_rna_qc.pct_counts_mt_max", float(pct_mt_upper),
         f"{pct_mt_k}*MAD above median(pct_counts_mt); clamped to [{pct_mt_floor}, {pct_mt_ceil}]"),
        ("s1_rna_qc.pct_counts_ribo_max", float(pct_ribo_max),
         "Soft ribosomal-protein ceiling (filters extreme stress/dying cells)."),
        ("s1_rna_qc.n_mt_genes_detected", n_mt,
         "Number of variables flagged as mitochondrial (sanity check for prefix detection)."),
        ("s1_rna_qc.n_ribo_genes_detected", n_ribo,
         "Number of variables flagged as ribosomal protein (Rps/Rpl/Mrps/Mrpl)."),
    ]:
        _prov.set_param(params_path, key, value,
                        source="derived", confidence="high", rationale=rat,
                        method={"name": "mad_thresholds" if "mt" in key or "min" in key or "max" in key else "var_prefix_detection",
                                "code_ref": "executor/methods/mad_thresholds.py"})

    # Apply filters
    keep = (
        (a.obs["total_counts"] >= c_lo)
        & (a.obs["total_counts"] <= c_hi)
        & (a.obs["n_genes_by_counts"] >= g_lo)
        & (a.obs["n_genes_by_counts"] <= g_hi)
        & (a.obs["pct_counts_mt"] <= pct_mt_upper)
        & (a.obs["pct_counts_ribo"] <= pct_ribo_max)
    )
    a_f = a[keep].copy()
    sc.pp.filter_genes(a_f, min_cells=min_cells)

    obs = a.obs
    cells_removed_per_metric = exclusive_removals({
        "total_counts": (
            (obs["total_counts"] >= c_lo) & (obs["total_counts"] <= c_hi)
        ).to_numpy(),
        "n_genes": (
            (obs["n_genes_by_counts"] >= g_lo) & (obs["n_genes_by_counts"] <= g_hi)
        ).to_numpy(),
        "pct_counts_mt": (obs["pct_counts_mt"] <= pct_mt_upper).to_numpy(),
        "pct_counts_ribo": (obs["pct_counts_ribo"] <= pct_ribo_max).to_numpy(),
    })

    # Save qc metrics pre/post
    _io.write_parquet_safe(a.obs, art / "qc_metrics_pre.parquet")
    _io.write_parquet_safe(a_f.obs, art / "qc_metrics_post.parquet")

    # QC violin figures are user-facing deliverables → checkpoint/qc_review/
    try:
        from .. import figures as _fig
        from ..run_paths import RunPaths
        figs_dir = RunPaths(run_dir).deliv_qc_review
        figs_dir.mkdir(parents=True, exist_ok=True)
        _fig.plot_qc_violin({
            "n_genes": a.obs["n_genes_by_counts"].to_numpy(),
            "total_counts": a.obs["total_counts"].to_numpy(),
            "pct_counts_mt": a.obs["pct_counts_mt"].to_numpy(),
            "pct_counts_ribo": a.obs["pct_counts_ribo"].to_numpy(),
        }, out_dir=figs_dir, stem="s1_rna_qc_violin_pre",
            title="RNA QC — pre-filter")
        _fig.plot_qc_violin({
            "n_genes": a_f.obs["n_genes_by_counts"].to_numpy(),
            "total_counts": a_f.obs["total_counts"].to_numpy(),
            "pct_counts_mt": a_f.obs["pct_counts_mt"].to_numpy(),
            "pct_counts_ribo": a_f.obs["pct_counts_ribo"].to_numpy(),
        }, out_dir=figs_dir, stem="s1_rna_qc_violin_post",
            title="RNA QC — post-filter")
    except Exception as e:
        log_event(run_dir, {"stage": "s1_rna_qc", "event": "plot_failed", "error": str(e)})

    result = {"n_cells_pre": int(a.n_obs), "n_cells_post": int(a_f.n_obs)}
    _io.write_text_safe(art / "qc_summary.json", json.dumps({
        "n_cells_pre": result["n_cells_pre"],
        "n_cells_post": result["n_cells_post"],
        "cells_removed_per_metric": cells_removed_per_metric,
    }, indent=2))
    _io.write_h5ad_safe(a_f, art / "rna_qc.h5ad")
    log_event(run_dir, {"stage": "s1_rna_qc", "event": "done",
                         "n_cells_pre": result["n_cells_pre"],
                         "n_cells_post": result["n_cells_post"],
                         "thresholds": {"total_counts": [c_lo, c_hi],
                                         "n_genes": [g_lo, g_hi],
                                         "pct_counts_mt_max": pct_mt_upper}})
    return result
