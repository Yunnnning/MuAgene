"""S1a — Ambient-RNA correction.

Whether correction runs is dataset- and plan-driven:

  - `method='none'` or nuclei sample type in the plan → pass-through.
  - `method='auto'` + raw input ref exists → SoupX (raw drops for soup profile).
  - `method='auto'` without raw matrix → DecontX (filtered counts only).
  - `method='decontx'` / `'soupx'` → forced backend (SoupX requires raw matrix).

R / DecontX / SoupX are required dependencies (`workflow/envs/processing.yaml`).
If they are missing when correction is requested, the stage raises and the run
fails — install or recreate the `muagene` env rather than silently skipping.
"""
from __future__ import annotations

import gc
import json
import shutil
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np

from .. import io as _io
from .. import provenance as _prov
from ..log import log_event
from ..methods import ambient as _amb


def _passthrough(src: Path, dst: Path) -> None:
    shutil.copy(src, dst)


def resolve_marker_genes(params_path: Path, plan_params: dict) -> list[str]:
    """Return marker gene list; parameters.yaml wins over frozen plan (like S2/S3)."""
    v = _prov.get_value(params_path, "s1a_ambient.marker_genes", None)
    if v:
        return list(v) if isinstance(v, list) else []
    entry = plan_params.get("marker_genes", {})
    mg = entry.get("value") if isinstance(entry, dict) else None
    return list(mg) if isinstance(mg, list) and mg else []


def _match_genes_to_var(adata: ad.AnnData, genes: list[str]) -> tuple[list[str], list[str]]:
    """Map requested symbols to var_names (case-insensitive); preserve request order."""
    lower_map = {str(v).lower(): str(v) for v in adata.var_names}
    found: list[str] = []
    missing: list[str] = []
    seen: set[str] = set()
    for g in genes:
        canon = lower_map.get(str(g).lower())
        if canon is None:
            missing.append(g)
        elif canon not in seen:
            found.append(canon)
            seen.add(canon)
    return found, missing


def _compute_tsne_from_layer(
    adata: ad.AnnData,
    layer: str,
    *,
    random_state: int = 42,
) -> np.ndarray:
    """t-SNE on one counts layer; used once and reused for before/after panels."""
    import scanpy as sc

    counts = adata.layers[layer]
    a = ad.AnnData(
        X=counts.copy() if hasattr(counts, "copy") else counts,
        obs=adata.obs,
        var=adata.var,
    )
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=min(2000, a.n_vars))
    a = a[:, a.var.highly_variable].copy()
    sc.pp.pca(a, n_comps=min(30, a.n_vars - 1))
    sc.tl.tsne(a, use_rep="X_pca", random_state=random_state, n_jobs=1)
    return np.asarray(a.obsm["X_tsne"])


def _log1p_expr_from_layer(
    adata: ad.AnnData,
    genes: list[str],
    layer: str,
    *,
    totals: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Log-normalized marker expression without copying the full matrix."""
    layer_mtx = adata.layers[layer]
    if totals is None:
        totals = np.asarray(layer_mtx.sum(axis=1)).ravel().astype(float)
    else:
        totals = np.asarray(totals, dtype=float).ravel()
    totals[totals == 0] = 1.0
    scale = 1e4 / totals
    found, _ = _match_genes_to_var(adata, genes)
    out: dict[str, np.ndarray] = {}
    for g in found:
        idx = adata.var_names.get_loc(g)
        x = layer_mtx[:, idx]
        raw = x.toarray().ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()
        out[g] = np.log1p(raw.astype(float) * scale)
    return out


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s1a_ambient"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"

    src = run_dir / "internal" / "artifacts" / "s0_ingest" / "rna_ingest.h5ad"
    s0_art = run_dir / "internal" / "artifacts" / "s0_ingest"
    dst = art / "rna_decontaminated.h5ad"

    # rna_ingest.h5ad is a deletable S0 *cache* — the post-QC cleanup removes it (S0's
    # DAG edge is its durable marker validation_report.json, which survives). It is a
    # deterministic function of the original input, so if it is absent (e.g. re-processing
    # a previously-approved run whose cache was cleaned) reconstruct it via the SSOT
    # io.load_rna_ingest and write it back — no full S0 re-run needed.
    if not src.exists():
        from ..run_paths import RunPaths
        import yaml as _yaml
        cfg = _yaml.safe_load(RunPaths(run_dir).run_yaml.read_text()) or {}
        rna_path = cfg.get("rna_path")
        if not rna_path:
            raise FileNotFoundError(
                f"s1a_ambient: S0 RNA ingest is missing at {src} and run.yaml has no "
                "`rna_path` to reconstruct it from. Re-run S0 ingest."
            )
        fmt = _prov.get_value(str(params_path), "ingest.rna_format", None)
        filtered_status = _prov.get_value(str(params_path), "ingest.rna_filtered_status", None)
        rebuilt, _raw, _diag = _io.load_rna_ingest(
            rna_path, fmt=fmt, filtered_status=filtered_status)
        src.parent.mkdir(parents=True, exist_ok=True)
        _io.write_h5ad_safe(rebuilt, src)
        log_event(run_dir, {"stage": "s1a_ambient", "event": "rna_ingest_reconstructed",
                            "rna_path": str(rna_path),
                            "note": ("rna_ingest.h5ad was absent (post-QC cleanup); rebuilt "
                                     "deterministically from the original input via "
                                     "io.load_rna_ingest — no S0 re-run")})

    p = plan["stages"].get("s1a_ambient", {}).get("parameters", {})
    method_choice = p.get("method", {}).get("value", "auto")
    max_contam = float(p.get("max_contamination", {}).get("value", 1.0))

    if method_choice == "none":
        chosen = "none"
    elif method_choice == "auto":
        chosen = "soupx" if _io.has_input_ref(s0_art, "rna_raw") else "decontx"
    elif method_choice in ("decontx", "soupx"):
        chosen = method_choice
    else:
        raise ValueError(
            f"s1a_ambient.method must be one of 'auto', 'none', 'decontx', "
            f"'soupx' — got {method_choice!r}."
        )

    if chosen == "none":
        _passthrough(src, dst)
        _prov.set_param(params_path, "s1a_ambient.method", "none",
                        source="user", confidence="high",
                        rationale="User-disabled ambient correction; S1a is a pass-through.")
        log_event(run_dir, {"stage": "s1a_ambient", "event": "passthrough",
                            "reason": "method=none"})
        n = _n_cells(dst)
        # summary.json is the durable stage-done marker (declared S1a output); it
        # must be written on every exit path, including pass-through, or status +
        # the S1a->S1 DAG edge break. rna_decontaminated.h5ad stays untracked.
        _io.write_text_safe(art / "summary.json", json.dumps(
            {"method": "none", "n_cells": n, "passthrough": True}, indent=2))
        return {"method": "none", "n_cells": n}

    a = ad.read_h5ad(src)
    if a.n_obs == 0 or a.n_vars == 0:
        _passthrough(src, dst)
        _prov.set_param(params_path, "s1a_ambient.method", "skipped_empty",
                        source="derived", confidence="high",
                        rationale="S0 RNA AnnData is empty (atac_only branch); pass-through.",
                        method={"name": "s1a.passthrough_empty",
                                "code_ref": "executor/stages/s1a_ambient.py"})
        _io.write_text_safe(art / "summary.json", json.dumps(
            {"method": "skipped_empty", "n_cells": 0, "passthrough": True}, indent=2))
        return {"method": "skipped_empty", "n_cells": 0}

    # Scratch dir for the SoupX/DecontX backend. It holds only transient working
    # files (no declared output, never read downstream) so it is removed as soon
    # as correction finishes — including on error — and a leftover from a prior
    # crashed run is cleared first so it cannot accumulate.
    work_dir = art / f"_work_{chosen}"
    if work_dir.exists():
        shutil.rmtree(work_dir, ignore_errors=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    try:
        if chosen == "soupx":
            if not _io.has_input_ref(s0_art, "rna_raw"):
                raise ValueError(
                    "SoupX requested but no raw RNA input ref exists — supply rna_raw_path "
                    "in run.yaml or set method=decontx/auto."
                )
            legacy_raw = s0_art / "rna_raw.h5ad"
            if legacy_raw.exists():
                raw_a = ad.read_h5ad(legacy_raw)
            else:
                raw_path, raw_fmt = _io.resolve_input_ref(s0_art / "rna_raw")
                raw_a = _io.load_rna(raw_path, fmt=raw_fmt)
            result = _amb.run_soupx(a, raw_a, work_dir=work_dir,
                                     max_contamination=max_contam)
            del raw_a
        else:
            result = _amb.run_decontx(a, work_dir=work_dir,
                                       max_contamination=max_contam)
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    a_corr = _amb.apply_correction(a, result)
    _wipe_marker_caches(art)
    _io.write_h5ad_safe(a_corr, dst)

    contam = np.asarray(result.contamination, dtype=float)
    contam_summary = {
        "method": result.method,
        "median_contamination": float(np.median(contam)) if contam.size else None,
        "mean_contamination": float(np.mean(contam)) if contam.size else None,
        "p90_contamination": float(np.quantile(contam, 0.90)) if contam.size else None,
        "n_high_contamination": int((contam > 0.20).sum()),
        "max_contamination_cap": max_contam,
    }
    _io.write_text_safe(art / "summary.json", json.dumps(
        {**contam_summary, **result.summary}, indent=2, default=str))
    import pandas as pd
    _io.write_parquet_safe(pd.DataFrame({"barcode": result.barcodes,
                  "contamination": contam}), art / "contamination.parquet")

    if method_choice in ("decontx", "soupx"):
        _prov_source = "user"
        _prov_rationale = (
            f"Explicit {result.method} choice from run.yaml/plan; "
            "per-cell rho written to obs['ambient_contamination'] and contamination.parquet."
        )
    else:  # method_choice == "auto"
        _prov_source = "derived"
        if chosen == "soupx":
            _prov_rationale = (
                f"Auto-selected {result.method} (raw RNA matrix present); "
                "per-cell rho written to obs['ambient_contamination'] and contamination.parquet."
            )
        else:
            _prov_rationale = (
                f"Auto-selected {result.method} (no raw matrix — filtered counts only); "
                "per-cell rho written to obs['ambient_contamination'] and contamination.parquet."
            )
    _prov.set_param(params_path, "s1a_ambient.method", result.method.lower(),
                    source=_prov_source, confidence="high",
                    rationale=_prov_rationale,
                    method={"name": f"ambient.run_{result.method.lower()}",
                            "code_ref": "executor/methods/ambient.py"})
    _prov.set_param(params_path, "s1a_ambient.median_contamination",
                    contam_summary["median_contamination"] or 0.0,
                    source="derived", confidence="high",
                    rationale=(f"Median rho across {len(contam)} cells "
                                f"({result.method})."),
                    method={"name": f"ambient.run_{result.method.lower()}",
                            "code_ref": "executor/methods/ambient.py"})
    _prov.set_param(params_path, "s1a_ambient.n_high_contamination_cells",
                    contam_summary["n_high_contamination"],
                    source="derived", confidence="high",
                    rationale="Cells with rho>0.20.",
                    method={"name": f"ambient.run_{result.method.lower()}",
                            "code_ref": "executor/methods/ambient.py"})

    method_label = result.method
    pre_counts = np.asarray(a_corr.layers["counts_raw"].sum(axis=1)).ravel()
    post_counts = np.asarray(a_corr.layers["counts"].sum(axis=1)).ravel()
    _write_cell_totals(art, a_corr.obs_names, pre_counts, post_counts)
    del a, result
    gc.collect()

    try:
        from .. import figures as _fig
        from ..run_paths import RunPaths
        figs_dir = RunPaths(run_dir).deliv_figures
        figs_dir.mkdir(parents=True, exist_ok=True)
        _fig.plot_counts_before_after(pre_counts, post_counts, out_dir=figs_dir,
                                       stem="s1a_ambient_counts_before_after",
                                       title=f"Total counts pre vs post {method_label}")
        marker_genes = resolve_marker_genes(params_path, p)
        if marker_genes:
            _plot_marker_genes(
                run_dir, marker_genes,
                force_tsne=False,
                write_params=False,
                refresh_qc=False,
                adata_in_memory=a_corr,
            )
    except Exception as e:
        log_event(run_dir, {"stage": "s1a_ambient", "event": "plot_failed", "error": str(e)})

    log_event(run_dir, {"stage": "s1a_ambient", "event": "done",
                        "method": method_label,
                        "median_contamination": contam_summary["median_contamination"],
                        "n_cells": int(a_corr.n_obs)})
    return {"method": method_label, "n_cells": int(a_corr.n_obs),
            "median_contamination": contam_summary["median_contamination"]}


_TSNE_CACHE_FILE = "tsne_coords_cache.parquet"
_CELL_TOTALS_FILE = "cell_totals.parquet"
_MARKER_CHECK_FILE = "marker_gene_check.json"


def _wipe_marker_caches(art_dir: Path) -> None:
    for name in (_TSNE_CACHE_FILE, _CELL_TOTALS_FILE, _MARKER_CHECK_FILE):
        (art_dir / name).unlink(missing_ok=True)


def _write_cell_totals(
    art_dir: Path,
    obs_names: list[str] | "pd.Index",
    totals_raw: np.ndarray,
    totals_corrected: np.ndarray,
) -> None:
    import pandas as pd

    pd.DataFrame({
        "obs_name": list(obs_names),
        "total_raw": np.asarray(totals_raw, dtype="float32"),
        "total_corrected": np.asarray(totals_corrected, dtype="float32"),
    }).to_parquet(art_dir / _CELL_TOTALS_FILE, index=False)


def _load_cell_totals(
    art_dir: Path,
    obs_names: list[str] | "pd.Index",
) -> tuple[np.ndarray, np.ndarray] | None:
    import pandas as pd

    path = art_dir / _CELL_TOTALS_FILE
    if not path.exists():
        return None
    try:
        cached = pd.read_parquet(path)
        if list(cached["obs_name"]) != list(obs_names):
            return None
        return (
            cached["total_raw"].to_numpy(dtype=float),
            cached["total_corrected"].to_numpy(dtype=float),
        )
    except Exception:
        return None


def _ensure_cell_totals(
    art_dir: Path,
    adata: ad.AnnData,
) -> tuple[np.ndarray, np.ndarray]:
    loaded = _load_cell_totals(art_dir, adata.obs_names)
    if loaded is not None:
        return loaded
    raw = np.asarray(adata.layers["counts_raw"].sum(axis=1)).ravel().astype(float)
    corr = np.asarray(adata.layers["counts"].sum(axis=1)).ravel().astype(float)
    _write_cell_totals(art_dir, adata.obs_names, raw, corr)
    return raw, corr


def _tsne_cache_matches(art_dir: Path, h5ad_path: Path) -> bool:
    cache = art_dir / _TSNE_CACHE_FILE
    if not cache.exists():
        return False
    try:
        import pandas as pd

        cached = pd.read_parquet(cache)
        a = ad.read_h5ad(h5ad_path, backed="r")
        try:
            return list(cached["obs_name"]) == list(a.obs_names)
        finally:
            if a.file is not None:
                a.file.close()
    except Exception:
        return False


def _load_cached_tsne(art_dir: Path) -> np.ndarray:
    import pandas as pd

    cached = pd.read_parquet(art_dir / _TSNE_CACHE_FILE)
    return cached[["tsne_x", "tsne_y"]].to_numpy()


def _load_or_compute_tsne(
    adata: ad.AnnData,
    art_dir: Path,
    *,
    random_state: int = 42,
    run_dir: Path | str | None = None,
    force: bool = False,
) -> np.ndarray:
    import pandas as pd

    cache = art_dir / _TSNE_CACHE_FILE
    if not force and cache.exists():
        try:
            cached = pd.read_parquet(cache)
            if list(cached["obs_name"]) == list(adata.obs_names):
                if run_dir is not None:
                    log_event(run_dir, {
                        "stage": "s1a_ambient", "event": "tsne_cache_hit",
                        "cache": str(cache),
                    })
                return cached[["tsne_x", "tsne_y"]].to_numpy()
        except Exception:
            pass

    reason = "force" if force else "miss"
    if run_dir is not None:
        log_event(run_dir, {
            "stage": "s1a_ambient", "event": "tsne_cache_miss",
            "reason": reason, "cache": str(cache),
        })

    work = adata.to_memory() if getattr(adata, "isbacked", False) else adata
    coords = _compute_tsne_from_layer(work, "counts_raw", random_state=random_state)
    pd.DataFrame({
        "obs_name": adata.obs_names,
        "tsne_x": coords[:, 0].astype("float32"),
        "tsne_y": coords[:, 1].astype("float32"),
    }).to_parquet(cache, index=False)
    return coords


def _write_marker_check_artifact(art: Path, found: list[str], missing: list[str]) -> None:
    _io.write_text_safe(
        art / _MARKER_CHECK_FILE,
        json.dumps({"found": found, "missing": missing}, indent=2),
    )


def _plot_marker_genes(
    run_dir: Path | str,
    marker_genes: list[str],
    *,
    force_tsne: bool = False,
    write_params: bool = True,
    refresh_qc: bool = False,
    adata_in_memory: ad.AnnData | None = None,
) -> dict[str, Any]:
    """Load data, extract expression, plot marker t-SNE figure; optionally refresh QC."""
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s1a_ambient"
    params_path = run_dir / "internal" / "parameters.yaml"
    dst = art / "rna_decontaminated.h5ad"

    if adata_in_memory is None and not dst.exists():
        raise FileNotFoundError(
            f"Post-correction RNA data not found: {dst}. "
            "S1a ambient correction must have completed before running the marker gene check."
        )

    owns_adata = adata_in_memory is None
    if adata_in_memory is not None:
        a_corr = adata_in_memory
        coords = _load_or_compute_tsne(a_corr, art, run_dir=run_dir, force=force_tsne)
    else:
        cache_hit = not force_tsne and _tsne_cache_matches(art, dst)
        if cache_hit:
            a_corr = ad.read_h5ad(dst, backed="r")
            log_event(run_dir, {
                "stage": "s1a_ambient", "event": "tsne_cache_hit",
                "cache": str(art / _TSNE_CACHE_FILE),
                "backed": True,
            })
            coords = _load_cached_tsne(art)
        else:
            a_corr = ad.read_h5ad(dst)
            coords = _load_or_compute_tsne(
                a_corr, art, run_dir=run_dir, force=force_tsne,
            )

    if "counts_raw" not in a_corr.layers or "counts" not in a_corr.layers:
        raise ValueError(
            "rna_decontaminated.h5ad is missing expected layers ('counts_raw', 'counts'). "
            "The file may have been produced by an older pipeline version."
        )

    found, missing = _match_genes_to_var(a_corr, marker_genes)
    totals_raw, totals_corr = _ensure_cell_totals(art, a_corr)
    expr_pre = _log1p_expr_from_layer(
        a_corr, marker_genes, "counts_raw", totals=totals_raw,
    )
    expr_post = _log1p_expr_from_layer(
        a_corr, marker_genes, "counts", totals=totals_corr,
    )

    if owns_adata and getattr(a_corr, "isbacked", False) and a_corr.file is not None:
        a_corr.file.close()
    if owns_adata:
        del a_corr
        gc.collect()

    if found:
        from .. import figures as _fig
        from ..run_paths import RunPaths
        figs_dir = RunPaths(run_dir).deliv_figures
        figs_dir.mkdir(parents=True, exist_ok=True)
        _fig.plot_marker_genes_tsne(
            coords, coords, expr_pre, expr_post,
            out_dir=figs_dir, stem="s1a_ambient_marker_genes",
            genes=found,
        )

    _write_marker_check_artifact(art, found, missing)

    if write_params:
        _prov.set_param(
            params_path, "s1a_ambient.marker_genes", marker_genes,
            source="user", confidence="high",
            rationale="Marker genes requested at QC review stage.",
        )
    elif not _prov.get_value(params_path, "s1a_ambient.marker_genes", None):
        _prov.set_param(
            params_path, "s1a_ambient.marker_genes", marker_genes,
            source="inferred", confidence="high",
            rationale="Marker genes from preprocessing plan; plotted during S1a.",
            method={"name": "s1a.resolve_marker_genes",
                    "code_ref": "executor/stages/s1a_ambient.py"},
        )
    _prov.set_param(
        params_path, "s1a_ambient.marker_genes_missing", missing,
        source="derived", confidence="high",
        rationale="Requested symbols absent from the expression matrix.",
        method={"name": "s1a.match_genes_to_var",
                "code_ref": "executor/stages/s1a_ambient.py"},
    )

    log_event(run_dir, {
        "stage": "s1a_ambient", "event": "marker_gene_check_done",
        "marker_genes": marker_genes, "found": found, "missing": missing,
    })

    if refresh_qc:
        from . import post_qc_review as _pqr
        _pqr.propose(run_dir)

    return {"found": found, "missing": missing}


def run_marker_gene_check(
    run_dir: Path | str,
    marker_genes: list[str],
    *,
    force_tsne: bool = False,
    refresh_qc: bool = True,
) -> dict[str, Any]:
    """Post-hoc marker gene check at QC review: plot and optionally refresh QC reports."""
    return _plot_marker_genes(
        run_dir, marker_genes,
        force_tsne=force_tsne,
        write_params=True,
        refresh_qc=refresh_qc,
    )


def _n_cells(p: Path) -> int:
    try:
        a = ad.read_h5ad(p, backed="r")
        try:
            return int(a.n_obs)
        finally:
            if a.file is not None:
                a.file.close()
    except Exception:
        return 0
