"""S6 — Dim reduction + neighbors (RNA PCA + neighbors; ATAC neighbors on LSI already present)."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import scanpy as sc

from .. import provenance as _prov
from ..log import log_event


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s6_dimred"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"
    branch = _prov.current_branch(str(params_path))

    has_rna = branch in ("paired", "separate", "rna_only")
    has_atac = branch in ("paired", "separate", "atac_only")

    s6_params = plan["stages"].get("s6_dimred", {}).get("parameters", {})
    n_pcs = int(s6_params.get("rna_n_pcs", {}).get("value", 30))
    n_neighbors = int(s6_params.get("n_neighbors", {}).get("value", 15))

    # --- RNA ---
    if has_rna:
        a = ad.read_h5ad(run_dir / "internal" / "artifacts" / "s4_rna_norm" / "rna_norm.h5ad")
        if "highly_variable" in a.var:
            sc.pp.pca(a, n_comps=n_pcs, use_highly_variable=True)
        else:
            sc.pp.pca(a, n_comps=n_pcs)
        sc.pp.neighbors(a, n_neighbors=n_neighbors, n_pcs=n_pcs)
        a.write_h5ad(art / "rna_dimred.h5ad")
        _prov.set_param(params_path, "s6_dimred.rna_n_pcs", n_pcs,
                        source="default", confidence="medium",
                        rationale="Plan default; elbow-based refinement skipped for MVP")
    else:
        # atac_only — produce an empty placeholder so the Snakemake output exists.
        import scipy.sparse as sp
        ad.AnnData(X=sp.csr_matrix((0, 0))).write_h5ad(art / "rna_dimred.h5ad")

    # --- ATAC ---
    if has_atac:
        import snapatac2 as snap
        atac_h5 = run_dir / "internal" / "artifacts" / "s5_atac_lsi" / "atac_lsi.h5ad"
        if not atac_h5.exists():
            atac_h5 = run_dir / "internal" / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"
        if atac_h5.exists():
            adata = snap.read(str(atac_h5))
            try:
                snap.pp.knn(adata, n_neighbors=n_neighbors)
            except Exception as e:
                log_event(run_dir, {"stage": "s6_dimred", "event": "atac_knn_failed", "error": str(e)})
            try:
                adata.close()
            except Exception:
                pass

    log_event(run_dir, {"stage": "s6_dimred", "event": "done",
                        "n_pcs": n_pcs, "n_neighbors": n_neighbors, "branch": branch})
    return {"n_pcs": n_pcs, "n_neighbors": n_neighbors, "branch": branch}
