"""S1a — Ambient-RNA correction.

Dispatch (driven by what S0 produced):

  - `rna_raw.h5ad` exists  → SoupX (uses raw drops to estimate the soup profile).
  - `rna_raw.h5ad` absent  → DecontX (filtered counts only).
  - method='none'          → pass-through (writes the S0 ingest unchanged so the
                              downstream DAG always has the same artifact).

R / DecontX / SoupX missing at runtime → falls back to pass-through with
`s1a_ambient.method = "skipped_no_r"` and a warning logged. The pipeline never
hard-fails on ambient correction; users without R still get a working run, with
provenance that explicitly documents the deviation.
"""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np

from .. import provenance as _prov
from ..log import log_event
from ..methods import ambient as _amb


def _passthrough(src: Path, dst: Path) -> None:
    shutil.copy(src, dst)


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s1a_ambient"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"

    src = run_dir / "internal" / "artifacts" / "s0_ingest" / "rna_ingest.h5ad"
    raw_src = run_dir / "internal" / "artifacts" / "s0_ingest" / "rna_raw.h5ad"
    dst = art / "rna_decontaminated.h5ad"

    p = plan["stages"].get("s1a_ambient", {}).get("parameters", {})
    method_choice = p.get("method", {}).get("value", "auto")
    max_contam = float(p.get("max_contamination", {}).get("value", 1.0))

    # --- Resolve method ---------------------------------------------------
    if method_choice == "none":
        chosen = "none"
    elif method_choice == "auto":
        chosen = "soupx" if raw_src.exists() else "decontx"
    elif method_choice in ("decontx", "soupx"):
        chosen = method_choice
    else:
        raise ValueError(
            f"s1a_ambient.method must be one of 'auto', 'none', 'decontx', "
            f"'soupx' — got {method_choice!r}."
        )

    # --- Pass-through path -----------------------------------------------
    if chosen == "none":
        _passthrough(src, dst)
        _prov.set_param(params_path, "s1a_ambient.method", "none",
                        source="user", confidence="high",
                        rationale="User-disabled ambient correction; S1a is a pass-through.")
        log_event(run_dir, {"stage": "s1a_ambient", "event": "passthrough",
                            "reason": "method=none"})
        return {"method": "none", "n_cells": _n_cells(dst)}

    # --- DecontX / SoupX -------------------------------------------------
    a = ad.read_h5ad(src)
    if a.n_obs == 0 or a.n_vars == 0:
        _passthrough(src, dst)
        _prov.set_param(params_path, "s1a_ambient.method", "skipped_empty",
                        source="derived", confidence="high",
                        rationale="S0 RNA AnnData is empty (atac_only branch); pass-through.",
                        method={"name": "s1a.passthrough_empty",
                                "code_ref": "executor/stages/s1a_ambient.py"})
        return {"method": "skipped_empty", "n_cells": 0}

    work_dir = art / f"_work_{chosen}"
    work_dir.mkdir(parents=True, exist_ok=True)

    try:
        if chosen == "soupx":
            if not raw_src.exists():
                raise _amb.AmbientUnavailable(
                    "SoupX requested but no rna_raw.h5ad exists — supply rna_raw_path "
                    "in run.yaml or set method=decontx/auto."
                )
            raw_a = ad.read_h5ad(raw_src)
            result = _amb.run_soupx(a, raw_a, work_dir=work_dir,
                                     max_contamination=max_contam)
        else:
            result = _amb.run_decontx(a, work_dir=work_dir,
                                       max_contamination=max_contam)
    except _amb.AmbientUnavailable as e:
        # Graceful degradation: never block the pipeline on a missing R dep.
        log_event(run_dir, {"stage": "s1a_ambient", "event": "skipped_no_r",
                            "method_attempted": chosen, "error": str(e)})
        _passthrough(src, dst)
        _prov.set_param(params_path, "s1a_ambient.method", "skipped_no_r",
                        source="derived", confidence="high",
                        rationale=(f"Attempted {chosen}; failed with: {str(e)[:300]}. "
                                    "S1a degraded to pass-through; downstream stages "
                                    "operate on uncorrected counts."),
                        method={"name": "s1a.skipped_no_r",
                                "code_ref": "executor/stages/s1a_ambient.py"})
        return {"method": "skipped_no_r", "error": str(e)}

    a_corr = _amb.apply_correction(a, result)
    a_corr.write_h5ad(dst)

    # --- Persist diagnostics + provenance --------------------------------
    contam = np.asarray(result.contamination, dtype=float)
    contam_summary = {
        "method": result.method,
        "median_contamination": float(np.median(contam)) if contam.size else None,
        "mean_contamination": float(np.mean(contam)) if contam.size else None,
        "p90_contamination": float(np.quantile(contam, 0.90)) if contam.size else None,
        "n_high_contamination": int((contam > 0.20).sum()),
        "max_contamination_cap": max_contam,
    }
    import json as _json
    (art / "summary.json").write_text(_json.dumps(
        {**contam_summary, **result.summary}, indent=2, default=str))
    import pandas as pd
    pd.DataFrame({"barcode": result.barcodes,
                  "contamination": contam}).to_parquet(art / "contamination.parquet")

    _prov.set_param(params_path, "s1a_ambient.method", result.method.lower(),
                    source="derived", confidence="high",
                    rationale=(f"Auto-selected {result.method} based on "
                                "presence of a raw matrix; per-cell rho written to "
                                "obs['ambient_contamination'] and contamination.parquet."),
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

    # --- Figures ---------------------------------------------------------
    try:
        from .. import figures as _fig
        from ..run_paths import RunPaths
        figs_dir = RunPaths(run_dir).deliv_figures
        figs_dir.mkdir(parents=True, exist_ok=True)
        _fig.plot_contamination_hist(contam, out_dir=figs_dir,
                                      stem="s1a_ambient_contamination_hist",
                                      title=f"Ambient RNA contamination ({result.method})")
        # Before/after counts comparison (per-cell totals).
        pre_counts = np.asarray(a.layers.get("counts", a.X).sum(axis=1)).ravel()
        post_counts = np.asarray(result.corrected_counts.sum(axis=1)).ravel()
        _fig.plot_counts_before_after(pre_counts, post_counts, out_dir=figs_dir,
                                       stem="s1a_ambient_counts_before_after",
                                       title=f"Total counts pre vs post {result.method}")
    except Exception as e:
        log_event(run_dir, {"stage": "s1a_ambient", "event": "plot_failed", "error": str(e)})

    log_event(run_dir, {"stage": "s1a_ambient", "event": "done",
                        "method": result.method,
                        "median_contamination": contam_summary["median_contamination"],
                        "n_cells": int(a_corr.n_obs)})
    return {"method": result.method, "n_cells": int(a_corr.n_obs),
            "median_contamination": contam_summary["median_contamination"]}


def _n_cells(p: Path) -> int:
    try:
        return int(ad.read_h5ad(p, backed="r").n_obs)
    except Exception:
        return 0
