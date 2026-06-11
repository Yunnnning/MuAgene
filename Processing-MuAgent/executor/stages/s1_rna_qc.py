"""S1 — RNA QC thresholds (MAD-based) + filtering."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import scanpy as sc

from ..methods import qc_thresholds as _qct
from ..methods.qc_filter_stats import marginal_removals
from ..methods.qc_metrics import compute_rna_qc_metrics
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

    # Mito + ribosomal flagging + per-cell QC metrics (shared with the pre-plan
    # QC exploration). Sanity-check by gene-count: a typical mammalian mt gene
    # set is ~13 genes; a typical ribosomal protein set is ~80–90 genes. If we
    # get zero mt genes the MT filter is effectively disabled — log a warning.
    qc_counts = compute_rna_qc_metrics(a)
    n_mt = qc_counts["n_mt"]
    n_ribo = qc_counts["n_ribo"]
    if n_mt == 0:
        log_event(run_dir, {"stage": "s1_rna_qc", "event": "no_mt_genes_detected",
                             "note": "var_names may be Ensembl IDs; pct_counts_mt will be 0 and the MT filter is effectively disabled"})

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

    # Derive thresholds (shared with the pre-plan QC exploration); absolute
    # floors win when MAD lower bounds fall below them.
    th = _qct.rna_thresholds(
        a.obs, k_mad=k_mad, pct_mt_k=pct_mt_k, pct_mt_ceiling=pct_mt_ceil,
        pct_mt_floor=pct_mt_floor, min_counts_floor=min_counts_floor,
        min_genes_floor=min_genes_floor,
    )
    c_lo, c_hi = th["total_counts_min"], th["total_counts_max"]
    g_lo, g_hi = th["n_genes_min"], th["n_genes_max"]
    pct_mt_upper = th["pct_counts_mt_max"]

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
    ]:
        _prov.set_param(params_path, key, value,
                        source="derived", confidence="high", rationale=rat,
                        method={"name": "mad_thresholds",
                                "code_ref": "executor/methods/mad_thresholds.py"})

    # Per-metric pass masks (shared with the pre-plan QC exploration).
    masks = _qct.rna_pass_masks(a.obs, th, pct_ribo_max=pct_ribo_max)

    # Apply filters
    keep = (
        masks["total_counts"]
        & masks["n_genes"]
        & masks["pct_counts_mt"]
        & masks["pct_counts_ribo"]
    )
    a_f = a[keep].copy()
    sc.pp.filter_genes(a_f, min_cells=min_cells)

    cells_removed_per_metric = marginal_removals(masks)

    # Save qc metrics pre/post
    _io.write_parquet_safe(a.obs, art / "qc_metrics_pre.parquet")
    _io.write_parquet_safe(a_f.obs, art / "qc_metrics_post.parquet")

    # QC violin figures are user-facing deliverables → deliverables/figures/
    try:
        from .. import figures as _fig
        from ..run_paths import RunPaths
        figs_dir = RunPaths(run_dir).deliv_figures
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
                         "n_mt_genes_detected": n_mt,
                         "n_ribo_genes_detected": n_ribo,
                         "thresholds": {"total_counts": [c_lo, c_hi],
                                         "n_genes": [g_lo, g_hi],
                                         "pct_counts_mt_max": pct_mt_upper}})
    return result
