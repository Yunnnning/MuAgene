"""S1 — RNA QC thresholds (MAD-based) + filtering."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import scanpy as sc

from ..methods import mad_thresholds as _mad
from .. import provenance as _prov
from ..log import log_event


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s1_rna_qc"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"

    in_path = run_dir / "internal" / "artifacts" / "s0_ingest" / "rna_ingest.h5ad"
    a = ad.read_h5ad(in_path)

    # Detect species-specific mito prefix
    mt_prefixes = ("mt-", "MT-", "Mt-")
    a.var["mt"] = a.var_names.str.startswith(tuple(mt_prefixes))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Parameters from plan
    params = plan["stages"]["s1_rna_qc"]["parameters"]
    k_mad = params["k_mad"]["value"]
    pct_mt_k = params["pct_mt_k"]["value"]
    pct_mt_ceil = params["pct_mt_ceiling"]["value"]
    pct_mt_floor = params["pct_mt_floor"]["value"]
    min_cells = int(params["min_cells_per_gene"]["value"])
    min_counts_floor = params["min_counts_floor"]["value"]

    # Apply floor before MAD
    keep_floor = a.obs["total_counts"] >= min_counts_floor
    a_for_mad = a[keep_floor].copy()

    # Derive thresholds
    c_lo, c_hi = _mad.log_mad_bounds(a_for_mad.obs["total_counts"].to_numpy(), k=k_mad)
    g_lo, g_hi = _mad.log_mad_bounds(a_for_mad.obs["n_genes_by_counts"].to_numpy(), k=k_mad)
    pct_mt_upper = _mad.upper_bound(a_for_mad.obs["pct_counts_mt"].to_numpy(),
                                     k=pct_mt_k, floor=pct_mt_floor, ceiling=pct_mt_ceil)

    # Record as provenance
    for key, value, rat in [
        ("s1_rna_qc.total_counts_min", float(c_lo), "MAD lower bound on log1p(total_counts)"),
        ("s1_rna_qc.total_counts_max", float(c_hi), "MAD upper bound on log1p(total_counts)"),
        ("s1_rna_qc.n_genes_min", float(g_lo), "MAD lower bound on log1p(n_genes)"),
        ("s1_rna_qc.n_genes_max", float(g_hi), "MAD upper bound on log1p(n_genes)"),
        ("s1_rna_qc.pct_counts_mt_max", float(pct_mt_upper),
         f"{pct_mt_k}*MAD above median(pct_counts_mt); clamped to [{pct_mt_floor}, {pct_mt_ceil}]"),
    ]:
        _prov.set_param(params_path, key, value,
                        source="derived", confidence="high", rationale=rat,
                        method={"name": "mad_thresholds",
                                "code_ref": "executor/methods/mad_thresholds.py"})

    # Apply filters
    keep = (
        (a.obs["total_counts"] >= c_lo)
        & (a.obs["total_counts"] <= c_hi)
        & (a.obs["n_genes_by_counts"] >= g_lo)
        & (a.obs["n_genes_by_counts"] <= g_hi)
        & (a.obs["pct_counts_mt"] <= pct_mt_upper)
    )
    a_f = a[keep].copy()
    sc.pp.filter_genes(a_f, min_cells=min_cells)

    # Save qc metrics pre/post
    a.obs.to_parquet(art / "qc_metrics_pre.parquet")
    a_f.obs.to_parquet(art / "qc_metrics_post.parquet")

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
        }, out_dir=figs_dir, stem="s1_rna_qc_violin_pre",
            title="RNA QC — pre-filter")
        _fig.plot_qc_violin({
            "n_genes": a_f.obs["n_genes_by_counts"].to_numpy(),
            "total_counts": a_f.obs["total_counts"].to_numpy(),
            "pct_counts_mt": a_f.obs["pct_counts_mt"].to_numpy(),
        }, out_dir=figs_dir, stem="s1_rna_qc_violin_post",
            title="RNA QC — post-filter")
    except Exception as e:
        log_event(run_dir, {"stage": "s1_rna_qc", "event": "plot_failed", "error": str(e)})

    a_f.write_h5ad(art / "rna_qc.h5ad")
    log_event(run_dir, {"stage": "s1_rna_qc", "event": "done",
                         "n_cells_pre": int(a.n_obs), "n_cells_post": int(a_f.n_obs),
                         "thresholds": {"total_counts": [c_lo, c_hi],
                                         "n_genes": [g_lo, g_hi],
                                         "pct_counts_mt_max": pct_mt_upper}})
    return {"n_cells_pre": int(a.n_obs), "n_cells_post": int(a_f.n_obs)}
