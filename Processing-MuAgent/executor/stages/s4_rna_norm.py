"""S4 — RNA normalization (log-normalize, target_sum=1e4) + HVG (seurat_v3 on counts)."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import scanpy as sc

from .. import provenance as _prov
from ..log import log_event


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s4_rna_norm"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"

    a = ad.read_h5ad(run_dir / "internal" / "artifacts" / "s3_doublets" / "rna_post_doublet.h5ad")

    p = plan["stages"]["s4_rna_norm"]["parameters"]
    target_sum = p["target_sum"]["value"]
    flavor = p["hvg_flavor"]["value"]
    n_top_cap = int(p["hvg_n_top_genes"]["value"])

    # Guard: seurat_v3 requires raw counts
    if flavor == "seurat_v3":
        if "counts" not in a.layers:
            raise RuntimeError("seurat_v3 HVG requires a 'counts' layer with raw integer counts")

    # HVG first (on counts), then normalize
    n_top = min(n_top_cap, max(200, int(0.1 * a.n_vars)))
    try:
        sc.pp.highly_variable_genes(a, n_top_genes=n_top, flavor=flavor, layer="counts")
    except Exception as e:
        log_event(run_dir, {"stage": "s4_rna_norm", "event": "hvg_seurat_v3_failed", "error": str(e)})
        # Fallback: normalize first, then cell-ranger flavor (less ideal)
        sc.pp.normalize_total(a, target_sum=target_sum)
        sc.pp.log1p(a)
        sc.pp.highly_variable_genes(a, n_top_genes=n_top, flavor="seurat")
    else:
        sc.pp.normalize_total(a, target_sum=target_sum)
        sc.pp.log1p(a)

    _prov.set_param(params_path, "s4_rna_norm.hvg_n_top_genes", int(n_top),
                    source="derived", confidence="high",
                    rationale=f"min({n_top_cap}, 0.1 * n_genes={int(0.1*a.n_vars)})",
                    method={"name": "hvg_cap", "code_ref": "executor/stages/s4_rna_norm.py"})

    a.write_h5ad(art / "rna_norm.h5ad")
    log_event(run_dir, {"stage": "s4_rna_norm", "event": "done",
                         "n_cells": int(a.n_obs), "n_hvg": int(a.var.get("highly_variable", 0).sum()) if "highly_variable" in a.var else None})
    return {"n_cells": int(a.n_obs)}
