"""S7 — Leiden clustering at fixed per-modality resolutions.

Clustering runs automatically with no resolution sweep and no user checkpoint:
the resolutions come from the plan parameters `s7_clustering.rna_resolution`
(default 0.7) and `s7_clustering.atac_resolution` (default 0.5). A user `revise`
recorded in parameters.yaml wins over the plan default (same overlay rule as the
QC stages). Leiden is run ONCE per modality and the applied resolutions are
recorded in parameters.yaml so the final review notebook and manifest can
surface them.

Outputs:
  - rna_clustered.h5ad      with final `leiden_rna` labels (empty stub when the
                            branch has no RNA modality)
  - atac_leiden_labels.parquet  barcode → `leiden_atac` (ATAC-bearing branches)
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from .. import io as _io
from .. import provenance as _prov
from ..log import log_event
from ..defaults import QC_DEFAULTS as _D


def _record_applied(params_path: Path, key: str, value: float) -> None:
    """Persist the applied resolution so the notebook/manifest can read it,
    without clobbering a user `revise` (which already wrote the key)."""
    if _prov.get_value(params_path, key, None) is None:
        _prov.set_param(params_path, key, value, source="default", confidence="high",
                        rationale="Fixed default Leiden resolution.")


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s7_clustering"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"
    branch = _prov.current_branch(str(params_path))
    has_rna = branch in ("paired", "separate", "rna_only")
    has_atac = branch in ("paired", "separate", "atac_only")

    plan_params = plan["stages"]["s7_clustering"]["parameters"]
    rna_res = float(_prov.effective_value(params_path, plan_params, "s7_clustering",
                                          "rna_resolution", _D["s7_clustering"]["rna_resolution"]))
    atac_res = float(_prov.effective_value(params_path, plan_params, "s7_clustering",
                                           "atac_resolution", _D["s7_clustering"]["atac_resolution"]))
    seed = int(_prov.effective_value(params_path, plan_params, "s7_clustering",
                                     "random_state", _D["s7_clustering"]["random_state"]))

    # --- RNA final labels ---
    if has_rna:
        _record_applied(params_path, "s7_clustering.rna_resolution", rna_res)
        rna = ad.read_h5ad(run_dir / "internal" / "artifacts" / "s6_neighbors" / "rna_neighbors.h5ad")
        sc.tl.leiden(rna, resolution=rna_res, random_state=seed, key_added="leiden_rna")
        _io.write_h5ad_safe(rna, art / "rna_clustered.h5ad")
    else:
        import scipy.sparse as sp
        _io.write_h5ad_safe(ad.AnnData(X=sp.csr_matrix((0, 0))), art / "rna_clustered.h5ad")

    # --- ATAC final labels ---
    if has_atac:
        _record_applied(params_path, "s7_clustering.atac_resolution", atac_res)
        try:
            import snapatac2 as snap
            atac_h5 = run_dir / "internal" / "artifacts" / "s5_atac_spectral" / "atac_spectral.h5ad"
            if not atac_h5.exists():
                atac_h5 = run_dir / "internal" / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"
            adata = snap.read(str(atac_h5))
            snap.tl.leiden(adata, resolution=atac_res, random_state=seed,
                           key_added="leiden_atac")
            try:
                leiden_col = adata.obs["leiden_atac"]
            except Exception:
                leiden_col = np.asarray(adata.obs["leiden_atac"])
            _io.write_parquet_safe(pd.DataFrame({
                "barcode": [str(x) for x in adata.obs_names],
                "leiden_atac": np.asarray(leiden_col).astype(str),
            }), art / "atac_leiden_labels.parquet", index=False)
            try:
                adata.close()
            except Exception:
                pass
        except Exception as e:
            log_event(run_dir, {"stage": "s7_clustering", "event": "atac_finalize_failed",
                                "error": str(e)})
            raise

    log_event(run_dir, {"stage": "s7_clustering", "event": "done",
                        "rna_resolution": rna_res if has_rna else None,
                        "atac_resolution": atac_res if has_atac else None})
    return {"rna_resolution": rna_res if has_rna else None,
            "atac_resolution": atac_res if has_atac else None}
