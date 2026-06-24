"""Pre-plan QC exploration.

Applies the **default** plan thresholds to the loaded data so the plan-review
document can show, per threshold, how many cells would be removed — plus
histograms with the cutoffs drawn on. This runs inside the merged S0 ingest job
(``executor.stages.s0_ingest``): the RNA matrix loaded there is passed in via
``rna_adata`` so it is never re-read, and the ATAC fragment import happens once
in that same job.

Per-cell QC metrics are persisted as small parquets (``rna_qc_metrics.parquet`` /
``atac_qc_metrics.parquet``) so a `revise` at the plan-review checkpoint can
re-derive thresholds, re-count removals, and re-draw histograms via
``rederive_from_metrics`` with no heavy reload.

Counting is non-exclusive: each threshold counts every cell that fails it,
evaluated independently on the full unfiltered dataset (see
``methods.qc_filter_stats.marginal_removals``). Threshold derivation is shared with
the real stages via ``methods.qc_thresholds`` so the preview matches S1/S2.

All figures are written to ``deliverables/figures/``; plan_review.md only
references them by relative path (the ``deliverables/plan/`` layout is unchanged).
FRiP is not computed here (needs MACS3); it is marked "computed at runtime".
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import numpy as np

from . import io as _io
from . import qc_tables as _qc_tables
from . import provenance as _prov
from .md_tables import fmt as _fmt
from .figures import (
    DEFAULT_MIN_COUNTS_FLOOR,
    DEFAULT_MIN_GENES_FLOOR,
    DEFAULT_N_FRAG_FLOOR,
    DEFAULT_N_FRAG_K_MAD,
    DEFAULT_NUC_MAX,
    DEFAULT_PCT_MT_CEILING,
    DEFAULT_PCT_MT_FLOOR,
    DEFAULT_PCT_MT_K,
    DEFAULT_PCT_MT_REFS,
    DEFAULT_PCT_RIBO_MAX,
    DEFAULT_TSS_MAX,
    DEFAULT_TSS_MIN,
    DEFAULT_TOTAL_COUNTS_K_MAD,
    DEFAULT_N_GENES_K_MAD,
    QC_EXPLORE_ATAC_TITLE,
    QC_EXPLORE_RNA_TITLE,
    RNA_HIST_EDGE_COLOR,
    RNA_HIST_FILL_ALPHA,
    RNA_HIST_FILL_COLOR,
    build_fixed_range_markers,
    build_mad_range_markers,
    build_upper_only_markers,
    default_atac_fragment_bounds,
    default_rna_thresholds,
    plot_qc_threshold_histograms,
    qc_hist_panel,
)
from .qc_tables import (
    SKIP_ABOVE_COUNTS,
    SKIP_ABOVE_GENES,
    SKIP_ABOVE_NUC,
    SKIP_ABOVE_TSS,
    SKIP_PCT_AT,
)
from .log import log_event
from .methods import qc_thresholds as _qct
from .methods.qc_filter_stats import marginal_removals
from .run_paths import RunPaths
from .defaults import QC_DEFAULTS as _D

S1_FIGURE_STEM = "s0_rna_data_explore"
S2_FIGURE_STEM = "s0_atac_data_explore"


def _pval(params: dict[str, Any], key: str, default: Any) -> Any:
    entry = params.get(key, {})
    if isinstance(entry, dict) and "value" in entry:
        return entry["value"]
    return default


def _effective_stage_params(run_dir: Path | str, plan: dict[str, Any], stage: str) -> dict[str, Any]:
    """Plan params for ``stage`` overlaid with parameters.yaml overrides.

    So the QC-exploration preview honours a ``revise`` at the plan-review
    checkpoint (the override wins over the frozen plan), matching what the real
    S1/S2 stages apply via ``provenance.effective_value``.
    """
    params_path = Path(run_dir) / "internal" / "parameters.yaml"
    plan_params = plan.get("stages", {}).get(stage, {}).get("parameters", {})
    return _prov.effective_params(params_path, plan_params, stage)


def _col_to_numpy(adata, key: str) -> np.ndarray:
    """SnapATAC2 obs is polars; convert one column to a float numpy array."""
    try:
        col = adata.obs[key]
    except Exception:
        return np.array([], dtype=float)
    try:
        arr = col.to_numpy()
    except AttributeError:
        arr = np.asarray(col)
    return np.asarray(arr, dtype=float)


# Per-cell QC metric columns persisted to parquet. Keeping these (rather than only
# the pass/fail counts) lets a `revise` at the plan-review checkpoint re-derive
# thresholds, re-count removals, and re-draw histograms with pure numpy — no h5ad
# reload, no fragment re-import. Histograms need the per-cell arrays anyway.
RNA_METRIC_COLS = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]
ATAC_METRIC_COLS = ["n_fragment", "tsse", "nucleosome_signal"]
RNA_METRICS_PARQUET = "rna_qc_metrics.parquet"
ATAC_METRICS_PARQUET = "atac_qc_metrics.parquet"

# The SnapATAC2-backed AnnData written by the explore import, plus a small meta
# sidecar describing it. s2_atac_qc reuses both (with the per-cell metrics
# parquet) so the heavy fragment import is not repeated.
ATAC_SNAP_EXPLORE_H5AD = "atac_snap_explore.h5ad"
ATAC_EXPLORE_META = "atac_explore_meta.json"
# Pre-filtering fragment-size distribution (SnapATAC2 frag_size_distr vector),
# persisted so the cheap re-derive path can redraw the QC grid's 4th panel
# without re-importing fragments.
ATAC_FRAG_SIZE_DISTR_NPY = "atac_frag_size_distr.npy"


# --- RNA -------------------------------------------------------------------

def _rna_qc_from_metrics(
    obs,
    params: dict[str, Any],
    figs_dir: Path,
) -> dict[str, Any]:
    """Derive RNA thresholds, count marginal removals, and draw histograms from a
    per-cell metrics frame (columns: RNA_METRIC_COLS). Pure compute + plotting —
    shared by the heavy explore path and the cheap re-derive path. ``params`` is
    the effective (override-overlaid) parameter set."""
    total_counts_k_mad = _pval(params, "total_counts_k_mad", DEFAULT_TOTAL_COUNTS_K_MAD)
    n_genes_k_mad = _pval(params, "n_genes_k_mad", DEFAULT_N_GENES_K_MAD)
    pct_mt_k = _pval(params, "pct_mt_k", DEFAULT_PCT_MT_K)
    pct_mt_ceil = _pval(params, "pct_mt_ceiling", DEFAULT_PCT_MT_CEILING)
    pct_mt_floor = _pval(params, "pct_mt_floor", DEFAULT_PCT_MT_FLOOR)
    min_counts_floor = _pval(params, "min_counts_floor", DEFAULT_MIN_COUNTS_FLOOR)
    min_genes_floor = float(_pval(params, "min_genes_floor", DEFAULT_MIN_GENES_FLOOR))
    pct_ribo_max = float(_pval(params, "pct_ribo_max", DEFAULT_PCT_RIBO_MAX))
    # Manual overrides pin the effective bound; the MAD/floor value is still
    # computed (th["*_derived"]) and shown grey under the red override.
    tc_min_ov = _pval(params, "total_counts_min_override", None)
    tc_max_ov = _pval(params, "total_counts_max_override", None)
    ng_min_ov = _pval(params, "n_genes_min_override", None)
    ng_max_ov = _pval(params, "n_genes_max_override", None)
    mt_max_ov = _pval(params, "pct_counts_mt_max_override", None)

    th = _qct.rna_thresholds(
        obs, total_counts_k_mad=total_counts_k_mad, n_genes_k_mad=n_genes_k_mad,
        pct_mt_k=pct_mt_k, pct_mt_ceiling=pct_mt_ceil,
        pct_mt_floor=pct_mt_floor, min_counts_floor=min_counts_floor,
        min_genes_floor=min_genes_floor,
        total_counts_min_override=tc_min_ov, total_counts_max_override=tc_max_ov,
        n_genes_min_override=ng_min_ov, n_genes_max_override=ng_max_ov,
        pct_counts_mt_max_override=mt_max_ov,
    )
    masks = _qct.rna_pass_masks(obs, th, pct_ribo_max=pct_ribo_max)
    cells_removed = marginal_removals(masks)
    th = {**th, "pct_counts_ribo_max": pct_ribo_max}

    default_th = default_rna_thresholds(obs)

    tc_markers, tc_lo, tc_hi = build_mad_range_markers(
        applied_lo=th["total_counts_min"],
        applied_hi=th["total_counts_max"],
        default_lo=default_th["total_counts_min"],
        default_hi=default_th["total_counts_max"],
        default_mad_lo_raw=default_th["total_counts_mad_lo_raw"],
        default_floor=DEFAULT_MIN_COUNTS_FLOOR,
        hi_skip_above=SKIP_ABOVE_COUNTS,
        log_axis=True,
        derived_lo=th["total_counts_min_derived"] if tc_min_ov is not None else None,
        derived_hi=th["total_counts_max_derived"] if tc_max_ov is not None else None,
    )
    ng_markers, ng_lo, ng_hi = build_mad_range_markers(
        applied_lo=th["n_genes_min"],
        applied_hi=th["n_genes_max"],
        default_lo=default_th["n_genes_min"],
        default_hi=default_th["n_genes_max"],
        default_mad_lo_raw=default_th["n_genes_mad_lo_raw"],
        default_floor=DEFAULT_MIN_GENES_FLOOR,
        hi_skip_above=SKIP_ABOVE_GENES,
        log_axis=True,
        derived_lo=th["n_genes_min_derived"] if ng_min_ov is not None else None,
        derived_hi=th["n_genes_max_derived"] if ng_max_ov is not None else None,
    )
    mt_markers, _, mt_hi = build_upper_only_markers(
        applied_hi=th["pct_counts_mt_max"],
        default_hi=default_th["pct_counts_mt_max"],
        default_mad_hi_raw=default_th["pct_counts_mt_mad_raw"],
        hi_skip_above=SKIP_PCT_AT,
        pct=True,
        default_fixed_refs=DEFAULT_PCT_MT_REFS,
        derived_hi=th["pct_counts_mt_max_derived"] if mt_max_ov is not None else None,
    )
    ribo_markers, _, ribo_hi = build_upper_only_markers(
        applied_hi=pct_ribo_max,
        default_hi=DEFAULT_PCT_RIBO_MAX,
        hi_skip_above=SKIP_PCT_AT,
        pct=True,
    )

    plot_qc_threshold_histograms(
        {
            "total_counts": qc_hist_panel(
                obs["total_counts"], tc_markers,
                filter_lo=tc_lo, filter_hi=tc_hi, log=True,
            ),
            "n_genes": qc_hist_panel(
                obs["n_genes_by_counts"], ng_markers,
                filter_lo=ng_lo, filter_hi=ng_hi, log=True,
            ),
            "pct_counts_mt": qc_hist_panel(
                obs["pct_counts_mt"], mt_markers, filter_hi=mt_hi,
            ),
            "pct_counts_ribo": qc_hist_panel(
                obs["pct_counts_ribo"], ribo_markers, filter_hi=ribo_hi,
            ),
        },
        out_dir=figs_dir, stem=S1_FIGURE_STEM,
        title=QC_EXPLORE_RNA_TITLE,
        fill_color=RNA_HIST_FILL_COLOR,
        edge_color=RNA_HIST_EDGE_COLOR,
        fill_alpha=RNA_HIST_FILL_ALPHA,
    )
    return {
        "thresholds": th,
        "cells_removed": cells_removed,
        "n_cells": int(len(obs["total_counts"])),
        "figure_stem": S1_FIGURE_STEM,
        "metrics_parquet": RNA_METRICS_PARQUET,
    }


def _explore_rna(run_dir: Path, plan: dict[str, Any], figs_dir: Path,
                 art: Path, adata=None) -> dict[str, Any] | None:
    """Compute per-cell RNA QC metrics, persist them, and produce the threshold
    preview. When ``adata`` is supplied (the merged S0 job passes its already-loaded
    matrix) the h5ad reload is skipped."""
    import pandas as pd

    if adata is None:
        rna_path = run_dir / "internal" / "artifacts" / "s0_ingest" / "rna_ingest.h5ad"
        if not rna_path.exists():
            return None
        import anndata as ad
        a = ad.read_h5ad(rna_path)
    else:
        a = adata
    if a.n_obs == 0 or a.n_vars == 0:
        return None

    from .methods.qc_metrics import compute_rna_qc_metrics

    # Mito/ribo flagging + per-cell QC metrics (shared with s1_rna_qc).
    compute_rna_qc_metrics(a)

    # Persist per-cell metrics for the cheap re-derive path.
    metrics = pd.DataFrame({c: np.asarray(a.obs[c], dtype=float) for c in RNA_METRIC_COLS})
    _io.write_parquet_safe(metrics, art / RNA_METRICS_PARQUET)

    params = _effective_stage_params(run_dir, plan, "s1_rna_qc")
    return _rna_qc_from_metrics(metrics, params, figs_dir)


# --- ATAC ------------------------------------------------------------------

def _atac_qc_from_metrics(
    n_frag_values: np.ndarray, tss_values: np.ndarray, ns_values: np.ndarray,
    n_pre: int, params: dict[str, Any], figs_dir: Path,
    frag_size_distr: "np.ndarray | None" = None,
) -> dict[str, Any]:
    """Derive ATAC thresholds, count marginal removals, and draw histograms from
    per-cell metric arrays. Shared by the heavy explore path and re-derive path.

    ``frag_size_distr`` (optional) is the pre-filtering SnapATAC2 fragment-size
    distribution; when present it fills the grid's 4th panel."""
    k_mad = _pval(params, "n_fragments_k_mad", DEFAULT_N_FRAG_K_MAD)
    n_frag_floor = _pval(params, "n_fragments_floor", DEFAULT_N_FRAG_FLOOR)
    tss_min = float(_pval(params, "tss_enrichment_min", DEFAULT_TSS_MIN))
    tss_max = float(_pval(params, "tss_enrichment_max", DEFAULT_TSS_MAX))
    nuc_max = float(_pval(params, "nucleosome_signal_max", DEFAULT_NUC_MAX))
    frip_min = float(_pval(params, "frip_min", _D["s2_atac_qc"]["frip_min"]))
    nf_min_ov = _pval(params, "n_fragments_min_override", None)
    nf_max_ov = _pval(params, "n_fragments_max_override", None)

    f_lo, f_hi, f_lo_mad_raw, (f_lo_derived, f_hi_derived) = _qct.atac_n_fragment_bounds(
        n_frag_values, k_mad=k_mad, n_frag_floor=n_frag_floor,
        n_fragments_min_override=nf_min_ov, n_fragments_max_override=nf_max_ov,
    )
    masks = _qct.atac_pass_masks(
        n_frag_values, tss_values, ns_values,
        f_lo=f_lo, f_hi=f_hi, tss_min=tss_min, tss_max=tss_max,
        nuc_max=nuc_max, n_pre=n_pre,
    )
    cells_removed = marginal_removals(masks)
    th = {
        "n_fragments_min": f_lo, "n_fragments_max": f_hi,
        "n_fragments_mad_lo_raw": f_lo_mad_raw,
        "tss_enrichment_min": tss_min, "tss_enrichment_max": tss_max,
        "nucleosome_signal_max": nuc_max, "frip_min": frip_min,
    }

    extra_panel = None
    if frag_size_distr is not None and np.asarray(frag_size_distr).size:
        extra_panel = {
            "distr": np.asarray(frag_size_distr, dtype=float),
            "title": "fragment size distribution (pre-filtering)",
        }
    default_f_lo, default_f_hi, default_f_mad = default_atac_fragment_bounds(n_frag_values)
    nf_markers, nf_lo, nf_hi = build_mad_range_markers(
        applied_lo=f_lo,
        applied_hi=f_hi,
        default_lo=default_f_lo,
        default_hi=default_f_hi,
        default_mad_lo_raw=default_f_mad,
        default_floor=DEFAULT_N_FRAG_FLOOR,
        hi_skip_above=SKIP_ABOVE_COUNTS,
        log_axis=True,
        derived_lo=f_lo_derived if nf_min_ov is not None else None,
        derived_hi=f_hi_derived if nf_max_ov is not None else None,
    )
    tss_markers, tss_lo, tss_hi = build_fixed_range_markers(
        applied_lo=tss_min,
        applied_hi=tss_max,
        default_lo=DEFAULT_TSS_MIN,
        default_hi=DEFAULT_TSS_MAX,
        hi_skip_above=SKIP_ABOVE_TSS,
    )
    nuc_markers, _, nuc_hi = build_upper_only_markers(
        applied_hi=nuc_max,
        default_hi=DEFAULT_NUC_MAX,
        hi_skip_above=SKIP_ABOVE_NUC,
    )

    plot_qc_threshold_histograms(
        {
            "n_fragments": qc_hist_panel(
                n_frag_values, nf_markers,
                filter_lo=nf_lo, filter_hi=nf_hi, log=True,
            ),
            "tss_enrichment": qc_hist_panel(
                tss_values, tss_markers,
                filter_lo=tss_lo, filter_hi=tss_hi,
            ),
            "nucleosome_signal": qc_hist_panel(
                ns_values, nuc_markers, filter_hi=nuc_hi,
            ),
        },
        out_dir=figs_dir, stem=S2_FIGURE_STEM,
        title=QC_EXPLORE_ATAC_TITLE,
        extra_panel=extra_panel,
        extra_panel_slot=1,
    )
    return {
        "thresholds": th,
        "cells_removed": cells_removed,
        "n_cells": n_pre,
        "figure_stem": S2_FIGURE_STEM,
        "frip_runtime": True,
        "metrics_parquet": ATAC_METRICS_PARQUET,
    }


def _explore_atac(run_dir: Path, plan: dict[str, Any], figs_dir: Path,
                  art: Path) -> dict[str, Any] | None:
    """Import ATAC fragments once, compute per-cell metrics, persist them, and
    produce the threshold preview. The heavy ``import_fragments`` now runs inside
    the merged S0 job (this function is called from S0's process)."""
    import pandas as pd

    meta_path = run_dir / "internal" / "artifacts" / "s0_ingest" / "atac_ingest.json"
    if not meta_path.exists():
        return None
    import snapatac2 as snap

    atac_meta = json.loads(meta_path.read_text())
    fragments_path = atac_meta["fragments_path"]
    params_path = run_dir / "internal" / "parameters.yaml"

    # Optional barcode-translation shim (paired branch) — mirrors S2 so the
    # whitelist (RNA-space barcodes) matches the imported fragments.
    translation_parquet = (run_dir / "internal" / "artifacts" / "s0_ingest"
                           / "barcode_translation.parquet")
    if translation_parquet.exists():
        from . import translation as _translation
        translated_path = art / "atac_fragments.translated.tsv.gz"
        if not translated_path.exists():
            table = _translation.load_translation_parquet(translation_parquet)
            _translation.translate_fragments_file(fragments_path, translated_path, table)
        fragments_path = str(translated_path)

    genome = _prov.get_value(str(params_path), "ingest.genome_assembly", None)
    genome_ref = getattr(snap.genome, genome, None) if genome else None
    if genome_ref is None:
        log_event(run_dir, {"stage": "qc_explore", "event": "atac_skipped_no_genome",
                            "genome": genome})
        return None

    # Normalise chrom naming (Ensembl "1"/"MT" → UCSC "chr1"/"chrM") + bounds-filter
    # via the shared helper — same code path as S2, so the two cannot drift. Robust
    # to a missing tabix/bgzip (Python peek for naming, gzip fallback for output);
    # raises on hard failure instead of silently importing un-renamed fragments,
    # which would match zero chromosomes and yield an empty ("(no data)") figure.
    cbf_path, add_chr_prefix = _io.prepare_fragments_for_snapatac(
        fragments_path, genome_ref, out_dir=art,
        log=lambda d: log_event(run_dir, {"stage": "qc_explore", **d}),
    )
    fragments_path = str(cbf_path)

    whitelist = atac_meta.get("cell_barcode_whitelist")
    h5_out = art / ATAC_SNAP_EXPLORE_H5AD
    Path(h5_out).unlink(missing_ok=True)  # idempotent re-runs
    snap_tmp = art / "snapatac2_tmp"
    snap_tmp.mkdir(exist_ok=True)
    adata = snap.pp.import_fragments(
        fragments_path, chrom_sizes=genome_ref, file=str(h5_out),
        sorted_by_barcode=False, whitelist=whitelist, tempdir=snap_tmp,
    )
    try:
        snap.metrics.tsse(adata, genome_ref)
    except Exception as e:
        log_event(run_dir, {"stage": "qc_explore", "event": "tsse_failed", "error": str(e)})

    # Dataset-level fragment-size distribution over ALL pre-filtering cells, for the
    # QC-exploration grid's 4th panel. Computed here (planning phase) and kept
    # independent of S2's post-QC distribution, which is over threshold-passing
    # cells only — the two are different processes and must not be coupled.
    frag_size_distr: np.ndarray | None = None
    try:
        snap.metrics.frag_size_distr(adata, max_recorded_size=1000)
        frag_size_distr = np.asarray(adata.uns["frag_size_distr"], dtype=float)
    except Exception as e:
        log_event(run_dir, {"stage": "qc_explore", "event": "frag_size_distr_failed",
                            "error": str(e)})

    n_frag_values = _col_to_numpy(adata, "n_fragment")
    tss_values = _col_to_numpy(adata, "tsse")
    cell_barcodes = list(adata.obs_names)
    n_pre = int(adata.n_obs)
    try:
        adata.close()
    except Exception:
        pass

    # Fail loud on a degenerate import: zero cells, or no cell with any fragment.
    # This almost always means the fragment chromosome names were not normalised to
    # the SnapATAC2 reference convention (Ensembl "1" vs UCSC "chr1") — an execution
    # error that must surface rather than be masked as an empty "(no data)" figure.
    _max_frag = float(np.max(n_frag_values)) if n_pre else 0.0
    if n_pre == 0 or _max_frag <= 0:
        raise RuntimeError(
            f"qc_explore ATAC import produced no usable cells (n_cells={n_pre}, "
            f"max n_fragment={_max_frag}). This usually means the fragment chromosome "
            f"names did not match the SnapATAC2 reference ({genome!r}) — check the "
            "Ensembl→UCSC chr-renaming and that bgzip/tabix (or the gzip fallback) ran. "
            "Refusing to emit an empty QC figure."
        )

    try:
        ns_values = _io.nucleosome_signal_per_cell(fragments_path, cell_barcodes)
    except Exception as e:
        log_event(run_dir, {"stage": "qc_explore", "event": "nuc_signal_failed",
                            "error": str(e), "fallback": "all_zeros"})
        ns_values = np.zeros(n_pre, dtype=float)

    # Persist per-cell metrics for the cheap re-derive path.
    metrics = pd.DataFrame({
        "n_fragment": np.asarray(n_frag_values, dtype=float),
        "tsse": np.asarray(tss_values, dtype=float),
        "nucleosome_signal": np.asarray(ns_values, dtype=float),
    })
    _io.write_parquet_safe(metrics, art / ATAC_METRICS_PARQUET)

    # Persist the pre-filtering fragment-size distribution for the cheap re-derive path.
    fsd_artifact: str | None = None
    if frag_size_distr is not None and frag_size_distr.size:
        np.save(art / ATAC_FRAG_SIZE_DISTR_NPY, frag_size_distr)
        fsd_artifact = ATAC_FRAG_SIZE_DISTR_NPY

    # Record the imported object + chrom-bound-filtered fragments so s2_atac_qc
    # can reuse them instead of re-importing the fragments.
    _io.write_text_safe(art / ATAC_EXPLORE_META, json.dumps({
        "atac_snap_h5ad": str(h5_out),
        "fragments_path": str(fragments_path),
        "add_chr_prefix": bool(add_chr_prefix),
        "genome": genome,
        "whitelist": whitelist,
        "n_cells": n_pre,
        "metrics_parquet": ATAC_METRICS_PARQUET,
        "frag_size_distr_npy": fsd_artifact,
    }, indent=2))

    params = _effective_stage_params(run_dir, plan, "s2_atac_qc")
    return _atac_qc_from_metrics(n_frag_values, tss_values, ns_values, n_pre,
                                 params, figs_dir, frag_size_distr=frag_size_distr)


def rederive_from_metrics(run_dir: Path | str) -> Path:
    """Cheap re-render: recompute thresholds + removal counts + histograms from the
    persisted per-cell metrics parquets, with no h5ad reload or fragment re-import.

    This is the path taken after a `revise` at the plan-review checkpoint and by the
    standalone `executor plan-review` re-render once the parquets exist.
    """
    import pandas as pd

    run_dir = Path(run_dir)
    paths = RunPaths(run_dir)
    art = paths.stage_dir("qc_explore")
    art.mkdir(parents=True, exist_ok=True)
    figs_dir = paths.deliv_figures
    figs_dir.mkdir(parents=True, exist_ok=True)

    plan_path = paths.preprocessing_plan
    plan = json.loads(plan_path.read_text()) if plan_path.exists() else {}
    stages = plan.get("stages", {})

    out: dict[str, Any] = {}
    rna_parquet = art / RNA_METRICS_PARQUET
    if "s1_rna_qc" in stages and rna_parquet.exists():
        try:
            metrics = pd.read_parquet(rna_parquet)
            params = _effective_stage_params(run_dir, plan, "s1_rna_qc")
            out["s1_rna_qc"] = _rna_qc_from_metrics(metrics, params, figs_dir)
        except Exception as e:
            log_event(run_dir, {"stage": "qc_explore", "event": "rna_rederive_failed",
                                "error": str(e)})

    atac_parquet = art / ATAC_METRICS_PARQUET
    if "s2_atac_qc" in stages and atac_parquet.exists():
        try:
            metrics = pd.read_parquet(atac_parquet)
            params = _effective_stage_params(run_dir, plan, "s2_atac_qc")
            # Reload the persisted pre-filtering fragment-size distribution so the
            # fragment-size panel re-renders without re-importing fragments.
            fsd_path = art / ATAC_FRAG_SIZE_DISTR_NPY
            frag_size_distr = np.load(fsd_path) if fsd_path.exists() else None
            out["s2_atac_qc"] = _atac_qc_from_metrics(
                metrics["n_fragment"].to_numpy(),
                metrics["tsse"].to_numpy(),
                metrics["nucleosome_signal"].to_numpy(),
                int(len(metrics)), params, figs_dir,
                frag_size_distr=frag_size_distr,
            )
        except Exception as e:
            log_event(run_dir, {"stage": "qc_explore", "event": "atac_rederive_failed",
                                "error": str(e)})

    out_path = art / "qc_explore.json"
    _io.write_text_safe(out_path, json.dumps(out, indent=2, default=str))
    log_event(run_dir, {"stage": "qc_explore", "event": "rederive_done",
                        "modalities": sorted(out.keys())})
    return out_path


def effective_thresholds(run_dir: Path | str, stage: str) -> dict[str, Any] | None:
    """Live QC thresholds for ``stage`` from the persisted per-cell metrics and the
    EFFECTIVE (override-overlaid) params — computed with the same ``_qct`` SSOT
    functions S1/S2 apply.

    Read-only and cheap: loads only the metrics parquet, draws no figures, writes no
    files. Used by the ``revise`` binding-constraint preview so it reflects the current
    ``parameters.yaml`` overlay rather than the frozen ``qc_explore.json`` snapshot
    (which is only refreshed at S0 / plan-review and therefore lags a post-QC revise).
    Returns ``None`` when the stage's metrics parquet is absent.
    """
    import pandas as pd

    run_dir = Path(run_dir)
    paths = RunPaths(run_dir)
    art = paths.stage_dir("qc_explore")
    plan_path = paths.preprocessing_plan
    plan = json.loads(plan_path.read_text()) if plan_path.exists() else {}

    if stage == "s1_rna_qc":
        pq = art / RNA_METRICS_PARQUET
        if not pq.exists():
            return None
        obs = pd.read_parquet(pq)
        params = _effective_stage_params(run_dir, plan, "s1_rna_qc")
        th = _qct.rna_thresholds(
            obs,
            total_counts_k_mad=_pval(params, "total_counts_k_mad", DEFAULT_TOTAL_COUNTS_K_MAD),
            n_genes_k_mad=_pval(params, "n_genes_k_mad", DEFAULT_N_GENES_K_MAD),
            pct_mt_k=_pval(params, "pct_mt_k", DEFAULT_PCT_MT_K),
            pct_mt_ceiling=_pval(params, "pct_mt_ceiling", DEFAULT_PCT_MT_CEILING),
            pct_mt_floor=_pval(params, "pct_mt_floor", DEFAULT_PCT_MT_FLOOR),
            min_counts_floor=_pval(params, "min_counts_floor", DEFAULT_MIN_COUNTS_FLOOR),
            min_genes_floor=float(_pval(params, "min_genes_floor", DEFAULT_MIN_GENES_FLOOR)),
            total_counts_min_override=_pval(params, "total_counts_min_override", None),
            total_counts_max_override=_pval(params, "total_counts_max_override", None),
            n_genes_min_override=_pval(params, "n_genes_min_override", None),
            n_genes_max_override=_pval(params, "n_genes_max_override", None),
            pct_counts_mt_max_override=_pval(params, "pct_counts_mt_max_override", None),
        )
        return {**th, "pct_counts_ribo_max": float(_pval(params, "pct_ribo_max", DEFAULT_PCT_RIBO_MAX))}

    if stage == "s2_atac_qc":
        pq = art / ATAC_METRICS_PARQUET
        if not pq.exists():
            return None
        m = pd.read_parquet(pq)
        params = _effective_stage_params(run_dir, plan, "s2_atac_qc")
        lo, hi, lo_raw, (lo_d, hi_d) = _qct.atac_n_fragment_bounds(
            m["n_fragment"].to_numpy(),
            k_mad=_pval(params, "n_fragments_k_mad", DEFAULT_N_FRAG_K_MAD),
            n_frag_floor=_pval(params, "n_fragments_floor", DEFAULT_N_FRAG_FLOOR),
            n_fragments_min_override=_pval(params, "n_fragments_min_override", None),
            n_fragments_max_override=_pval(params, "n_fragments_max_override", None),
        )
        return {
            "n_fragments_min": lo, "n_fragments_max": hi, "n_fragments_mad_lo_raw": lo_raw,
            "tss_enrichment_min": float(_pval(params, "tss_enrichment_min", DEFAULT_TSS_MIN)),
            "tss_enrichment_max": float(_pval(params, "tss_enrichment_max", DEFAULT_TSS_MAX)),
            "nucleosome_signal_max": float(_pval(params, "nucleosome_signal_max", DEFAULT_NUC_MAX)),
            "frip_min": float(_pval(params, "frip_min", _D["s2_atac_qc"]["frip_min"])),
        }
    return None


def run(run_dir: Path | str, *, rna_adata=None) -> Path:
    """Run QC exploration. The merged S0 job passes its already-loaded RNA matrix
    via ``rna_adata`` so no reload happens. When called standalone (e.g. the
    plan-review re-render) and the per-cell metrics parquets already exist, the
    cheap re-derive path is used instead of any heavy load.
    """
    run_dir = Path(run_dir)
    paths = RunPaths(run_dir)
    art = paths.stage_dir("qc_explore")
    art.mkdir(parents=True, exist_ok=True)
    figs_dir = paths.deliv_figures
    figs_dir.mkdir(parents=True, exist_ok=True)

    plan_path = paths.preprocessing_plan
    plan = json.loads(plan_path.read_text()) if plan_path.exists() else {}
    stages = plan.get("stages", {})

    # Cheap path: no fresh in-memory matrix, and every per-cell metrics parquet the
    # active branch needs already exists — re-derive without any heavy load.
    if rna_adata is None:
        need_rna = "s1_rna_qc" in stages
        need_atac = "s2_atac_qc" in stages
        have_rna = (art / RNA_METRICS_PARQUET).exists()
        have_atac = (art / ATAC_METRICS_PARQUET).exists()
        if (have_rna or have_atac) and (not need_rna or have_rna) and (not need_atac or have_atac):
            return rederive_from_metrics(run_dir)

    out: dict[str, Any] = {}
    if "s1_rna_qc" in stages:
        try:
            rna = _explore_rna(run_dir, plan, figs_dir, art, adata=rna_adata)
            if rna is not None:
                out["s1_rna_qc"] = rna
        except Exception as e:
            log_event(run_dir, {"stage": "qc_explore", "event": "rna_failed", "error": str(e)})
    if "s2_atac_qc" in stages:
        try:
            atac = _explore_atac(run_dir, plan, figs_dir, art)
            if atac is not None:
                out["s2_atac_qc"] = atac
        except Exception as e:
            log_event(run_dir, {"stage": "qc_explore", "event": "atac_failed", "error": str(e)})

    out_path = art / "qc_explore.json"
    _io.write_text_safe(out_path, json.dumps(out, indent=2, default=str))
    log_event(run_dir, {"stage": "qc_explore", "event": "done",
                        "modalities": sorted(out.keys())})
    return out_path


# --- Rendering -------------------------------------------------------------

def _rna_table(data: dict[str, Any]) -> str:
    return _qc_tables.rna_removal_table(
        data["thresholds"], data["cells_removed"],
        include_note=True,
    )


def _atac_table(data: dict[str, Any]) -> str:
    return _qc_tables.atac_removal_table(
        data["thresholds"], data["cells_removed"],
        include_note=True,
        frip_runtime_note=True,
        frip_removed="—",
    )


def render_appendix_blocks(run_dir: Path | str) -> dict[str, str]:
    """Return ``{stage_id: markdown}`` (table + embedded histogram) for the
    s1_rna_qc / s2_atac_qc appendix sections, or ``{}`` when exploration is
    absent. Image src is relative to plan_review.md's directory."""
    paths = RunPaths(Path(run_dir))
    json_path = paths.artifact("qc_explore", "qc_explore.json")
    if not json_path.exists():
        return {}
    try:
        data = json.loads(json_path.read_text())
    except Exception:
        return {}

    md_parent = paths.plan_review_md.parent
    blocks: dict[str, str] = {}
    builders = {"s1_rna_qc": _rna_table, "s2_atac_qc": _atac_table}
    for stage, builder in builders.items():
        d = data.get(stage)
        if not d:
            continue
        n_cells = d.get("n_cells")
        intro = (
            f"Default thresholds applied to the {n_cells:,} loaded cells; each row "
            "counts every cell failing that threshold (independently)."
            if isinstance(n_cells, int) else
            "Default thresholds applied to the loaded data; each row counts every "
            "cell failing that threshold (independently)."
        )
        parts = [intro, "", builder(d)]
        stem = d.get("figure_stem")
        if stem:
            png = paths.deliv_figures_path(stem)
            if png.is_file():
                src = Path(os.path.relpath(png, md_parent)).as_posix()
                parts += ["", f"![{stage} QC exploration]({src})"]
        blocks[stage] = "\n".join(parts)
    return blocks
