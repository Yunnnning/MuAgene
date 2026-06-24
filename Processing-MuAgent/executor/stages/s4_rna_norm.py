"""S4 — RNA normalization (log-normalize, target_sum=1e4) + HVG (seurat_v3 on counts)."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import scanpy as sc

from .. import io as _io
from .. import provenance as _prov
from ..log import log_event


def _load_rna_postqc(run_dir: Path):
    """Read the post-QC, post-doublet RNA AnnData (carries ``layers['counts']``).

    Canonical source: the ``rna`` modality of the post-QC handoff h5mu
    (``deliverables/qc/post_qc_<run>.h5mu``, written by qc_handoff) — read with
    ``mudata.read_h5ad(path, "rna")`` so the atac modality is not loaded. Legacy /
    transition fallback: the transient ``s3_doublets/rna_post_doublet.h5ad`` when the
    h5mu is absent (a pre-dedup run, or qc_handoff not yet run).
    """
    from ..run_paths import RunPaths
    h5mu = RunPaths(run_dir).post_qc_h5mu
    if h5mu.exists():
        import mudata as mu
        return mu.read_h5ad(str(h5mu), "rna")
    legacy = run_dir / "internal" / "artifacts" / "s3_doublets" / "rna_post_doublet.h5ad"
    if legacy.exists():
        return ad.read_h5ad(legacy)
    raise FileNotFoundError(
        f"s4_rna_norm: post-QC RNA not found — neither {h5mu} nor {legacy} exists. "
        "Run qc_handoff (after post_qc_review approval) before S4."
    )


def _write_marker(art: Path, payload: dict[str, Any]) -> None:
    """Write the durable stage-done marker (norm_summary.json).

    norm_summary.json is the SOLE Snakemake-declared output and the durable
    stage-done marker (status + the S4 -> S6 dependency edge key off it). The large
    rna_norm.h5ad is written as an UNTRACKED working file: it is read by path by S6
    and removed by `finish-cleanup` once the run's processed deliverable exists, so
    keeping it out of the declared DAG means deleting it never triggers a re-run.
    Must be the LAST write so marker-exists <=> stage-done.
    """
    import json
    _io.write_text_safe(art / "norm_summary.json",
                        json.dumps({"stage": "s4_rna_norm", **payload}, indent=2))


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s4_rna_norm"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"
    branch = _prov.current_branch(str(params_path))

    if branch == "atac_only":
        import scipy.sparse as sp
        _io.write_h5ad_safe(ad.AnnData(X=sp.csr_matrix((0, 0))), art / "rna_norm.h5ad")
        log_event(run_dir, {"stage": "s4_rna_norm", "event": "skipped_no_rna",
                            "branch": branch})
        _write_marker(art, {"branch": branch, "skipped": True, "n_cells": 0})
        return {"n_cells": 0, "branch": branch}

    a = _load_rna_postqc(run_dir)

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

    _io.write_h5ad_safe(a, art / "rna_norm.h5ad")
    n_hvg = int(a.var["highly_variable"].sum()) if "highly_variable" in a.var else None
    log_event(run_dir, {"stage": "s4_rna_norm", "event": "done",
                         "n_cells": int(a.n_obs), "n_hvg": n_hvg})
    _write_marker(art, {"branch": branch, "n_cells": int(a.n_obs), "n_hvg": n_hvg})
    return {"n_cells": int(a.n_obs)}
