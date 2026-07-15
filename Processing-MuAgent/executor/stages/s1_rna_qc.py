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
from ..defaults import QC_DEFAULTS as _D


def _resolve_param(params_path: Path, plan_params: dict, name: str, default: Any = None) -> Any:
    """parameters.yaml wins over plan (so `executor revise` takes effect on re-run)."""
    return _prov.effective_value(params_path, plan_params, "s1_rna_qc", name, default)


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
    if branch in ("paired", "unpaired", "rna_only") and not s1a_path.exists():
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

    # Parameters from plan, overlaid with any parameters.yaml override (a user
    # `revise` of a recipe knob wins over the frozen plan — same rule as S2/S3).
    params = plan["stages"]["s1_rna_qc"]["parameters"]
    total_counts_k_mad = _resolve_param(params_path, params, "total_counts_k_mad", _D["s1_rna_qc"]["total_counts_k_mad"])
    n_genes_k_mad = _resolve_param(params_path, params, "n_genes_k_mad", _D["s1_rna_qc"]["n_genes_k_mad"])
    pct_mt_k = _resolve_param(params_path, params, "pct_mt_k", _D["s1_rna_qc"]["pct_mt_k"])
    pct_mt_ceil = _resolve_param(params_path, params, "pct_mt_ceiling", _D["s1_rna_qc"]["pct_mt_ceiling"])
    pct_mt_floor = _resolve_param(params_path, params, "pct_mt_floor", _D["s1_rna_qc"]["pct_mt_floor"])
    min_cells = int(_resolve_param(params_path, params, "min_cells_per_gene", _D["s1_rna_qc"]["min_cells_per_gene"]))
    min_counts_floor = _resolve_param(params_path, params, "min_counts_floor", _D["s1_rna_qc"]["min_counts_floor"])
    min_genes_floor = float(_resolve_param(params_path, params, "min_genes_floor", _D["s1_rna_qc"]["min_genes_floor"]))
    # Ribosomal upper-bound is recommended (not enforced strictly): some
    # tissues legitimately have very high ribo-protein expression.
    pct_ribo_max = float(_resolve_param(params_path, params, "pct_ribo_max", _D["s1_rna_qc"]["pct_ribo_max"]))
    # Manual overrides pin the effective MAD-derived bound to an exact value (the
    # MAD/floor derivation still runs and is recorded as the rationale + grey
    # reference line). Absent → derived behaviour unchanged.
    tc_min_ov = _resolve_param(params_path, params, "total_counts_min_override", None)
    tc_max_ov = _resolve_param(params_path, params, "total_counts_max_override", None)
    ng_min_ov = _resolve_param(params_path, params, "n_genes_min_override", None)
    ng_max_ov = _resolve_param(params_path, params, "n_genes_max_override", None)
    mt_max_ov = _resolve_param(params_path, params, "pct_counts_mt_max_override", None)

    # Derive thresholds (shared with the pre-plan QC exploration); absolute
    # floors win when MAD lower bounds fall below them.
    th = _qct.rna_thresholds(
        a.obs, total_counts_k_mad=total_counts_k_mad, n_genes_k_mad=n_genes_k_mad,
        pct_mt_k=pct_mt_k, pct_mt_ceiling=pct_mt_ceil,
        pct_mt_floor=pct_mt_floor, min_counts_floor=min_counts_floor,
        min_genes_floor=min_genes_floor,
        total_counts_min_override=tc_min_ov, total_counts_max_override=tc_max_ov,
        n_genes_min_override=ng_min_ov, n_genes_max_override=ng_max_ov,
        pct_counts_mt_max_override=mt_max_ov,
    )
    c_lo, c_hi = th["total_counts_min"], th["total_counts_max"]
    g_lo, g_hi = th["n_genes_min"], th["n_genes_max"]
    pct_mt_upper = th["pct_counts_mt_max"]

    # Warn when an override is more permissive than its recommended fixed bound
    # (override still wins — the user is trusted — but the deviation is surfaced
    # in the QC review report and the log).
    override_warnings: list[str] = []
    if tc_min_ov is not None and float(tc_min_ov) < float(min_counts_floor):
        override_warnings.append(
            f"total_counts lower bound override {float(tc_min_ov):.4g} is below the "
            f"recommended floor min_counts_floor={float(min_counts_floor):.4g}")
    if ng_min_ov is not None and float(ng_min_ov) < float(min_genes_floor):
        override_warnings.append(
            f"n_genes lower bound override {float(ng_min_ov):.4g} is below the "
            f"recommended floor min_genes_floor={float(min_genes_floor):.4g}")
    if mt_max_ov is not None and float(mt_max_ov) > float(pct_mt_ceil):
        override_warnings.append(
            f"pct_counts_mt upper bound override {float(mt_max_ov):.4g}% is above the "
            f"recommended ceiling pct_mt_ceiling={float(pct_mt_ceil):.4g}%")
    if override_warnings:
        log_event(run_dir, {"stage": "s1_rna_qc", "event": "override_below_floor",
                            "warnings": override_warnings})

    # Record as provenance. An active override is the user's choice (source=user,
    # no method); otherwise the bound is MAD/floor-derived (source=derived, method).
    _MAD_METHOD = {"name": "mad_thresholds", "code_ref": "executor/methods/mad_thresholds.py"}
    for key, value, override, derived, rat in [
        ("s1_rna_qc.total_counts_min", float(c_lo), tc_min_ov, th["total_counts_min_derived"],
         f"max(MAD lower bound on log1p(total_counts), min_counts_floor={min_counts_floor})"),
        ("s1_rna_qc.total_counts_max", float(c_hi), tc_max_ov, th["total_counts_max_derived"],
         f"MAD upper bound on log1p(total_counts) (total_counts_k_mad={total_counts_k_mad})"),
        ("s1_rna_qc.n_genes_min", float(g_lo), ng_min_ov, th["n_genes_min_derived"],
         f"max(MAD lower bound on log1p(n_genes), min_genes_floor={min_genes_floor})"),
        ("s1_rna_qc.n_genes_max", float(g_hi), ng_max_ov, th["n_genes_max_derived"],
         f"MAD upper bound on log1p(n_genes) (n_genes_k_mad={n_genes_k_mad})"),
        ("s1_rna_qc.pct_counts_mt_max", float(pct_mt_upper), mt_max_ov, th["pct_counts_mt_max_derived"],
         f"{pct_mt_k}*MAD above median(pct_counts_mt); clamped to [{pct_mt_floor}, {pct_mt_ceil}]"),
    ]:
        if override is not None:
            _prov.set_param(params_path, key, value, source="user", confidence="high",
                            rationale=f"Manual override (was MAD-derived {float(derived):.4g})")
        else:
            _prov.set_param(params_path, key, value, source="derived", confidence="high",
                            rationale=rat, method=_MAD_METHOD)
    _prov.set_param(params_path, "s1_rna_qc.pct_counts_ribo_max", float(pct_ribo_max),
                    source="derived", confidence="high",
                    rationale="Soft ribosomal-protein ceiling (filters extreme stress/dying cells).",
                    method=_MAD_METHOD)

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
            title="RNA violin plots (pre-filtering)")
        _fig.plot_qc_violin({
            "n_genes": a_f.obs["n_genes_by_counts"].to_numpy(),
            "total_counts": a_f.obs["total_counts"].to_numpy(),
            "pct_counts_mt": a_f.obs["pct_counts_mt"].to_numpy(),
            "pct_counts_ribo": a_f.obs["pct_counts_ribo"].to_numpy(),
        }, out_dir=figs_dir, stem="s1_rna_qc_violin_post",
            title="RNA violin plots (post-filtering)")
    except Exception as e:
        log_event(run_dir, {"stage": "s1_rna_qc", "event": "plot_failed", "error": str(e)})

    result = {"n_cells_pre": int(a.n_obs), "n_cells_post": int(a_f.n_obs)}
    _io.write_text_safe(art / "qc_summary.json", json.dumps({
        "n_cells_pre": result["n_cells_pre"],
        "n_cells_post": result["n_cells_post"],
        "cells_removed_per_metric": cells_removed_per_metric,
        "override_warnings": override_warnings,
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
