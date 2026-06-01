"""S6 — Dim reduction + neighbors (RNA PCA + neighbors; ATAC neighbors on spectral embedding from S5).

Standard scanpy preprocessing path (rev. 2026-04):
    log-normalize → HVG → optional sc.pp.scale → PCA → neighbors

`rna_scale` (plan param, default True) toggles the scaling step. `rna_n_pcs`
defaults to `"auto"` and is resolved via a chord-distance knee on the
cumulative explained-variance curve, capped at `rna_n_pcs_max`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import scanpy as sc

from .. import io as _io
from .. import provenance as _prov
from ..log import log_event


def _pick_pca_elbow(variance_ratio: np.ndarray, *, max_n: int) -> tuple[int, str]:
    """Return the elbow `n_pcs` from the cumulative explained-variance curve.

    Uses the chord-distance heuristic: the rank that is farthest from the
    chord between (rank=1, cumvar[0]) and (rank=N, cumvar[-1]) is the knee.
    Falls back to `min(max_n, len)` if the curve is too flat to admit a knee.
    """
    var = np.asarray(variance_ratio, dtype=float)
    var = var[np.isfinite(var)]
    if var.size <= 2:
        return int(min(max_n, max(var.size, 1))), "fallback (too few PCs)"
    cum = np.cumsum(var)
    n = cum.size
    x = np.arange(1, n + 1, dtype=float)
    p1 = np.array([x[0], cum[0]])
    p2 = np.array([x[-1], cum[-1]])
    chord = p2 - p1
    chord_norm = float(np.linalg.norm(chord))
    if chord_norm == 0.0:
        return int(min(max_n, n)), "fallback (degenerate chord)"
    chord_unit = chord / chord_norm
    rel = np.column_stack([x, cum]) - p1
    proj = rel - np.outer(rel @ chord_unit, chord_unit)
    distances = np.linalg.norm(proj, axis=1)
    knee = int(np.argmax(distances)) + 1
    knee = max(2, min(knee, max_n, n))
    return knee, f"chord-distance knee at PC{knee} (max search n={n})"


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s6_dimred"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"
    branch = _prov.current_branch(str(params_path))

    has_rna = branch in ("paired", "separate", "rna_only")
    has_atac = branch in ("paired", "separate", "atac_only")

    s6_params = plan["stages"].get("s6_dimred", {}).get("parameters", {})
    n_pcs_param = s6_params.get("rna_n_pcs", {}).get("value", "auto")
    n_pcs_max = int(s6_params.get("rna_n_pcs_max", {}).get("value", 50))
    do_scale = bool(s6_params.get("rna_scale", {}).get("value", True))
    n_neighbors = int(s6_params.get("n_neighbors", {}).get("value", 15))

    # --- RNA ---
    if has_rna:
        a = ad.read_h5ad(run_dir / "internal" / "artifacts" / "s4_rna_norm" / "rna_norm.h5ad")

        # Subset to HVGs before scaling/PCA to avoid densifying the full matrix.
        # sc.pp.scale with zero_center=True (default) creates a dense copy;
        # operating on the HVG subset (typically 2K genes) keeps peak memory low.
        if "highly_variable" in a.var:
            a_pca = a[:, a.var["highly_variable"]].copy()
        else:
            a_pca = a

        if do_scale:
            try:
                sc.pp.scale(a_pca, max_value=10)
            except Exception as e:
                log_event(run_dir, {"stage": "s6_dimred", "event": "scale_failed",
                                    "error": str(e)})

        # Elbow-resolved n_pcs. Compute up to n_pcs_max then trim to the elbow.
        n_seed = int(min(n_pcs_max, max(2, a_pca.n_vars - 1)))
        sc.pp.pca(a_pca, n_comps=n_seed)

        # Copy PCA results back to the full object (neighbors uses obsm["X_pca"]).
        a.obsm["X_pca"] = a_pca.obsm["X_pca"]
        a.uns["pca"] = a_pca.uns.get("pca", {})

        if isinstance(n_pcs_param, str) and n_pcs_param.strip().lower() == "auto":
            vr = np.asarray(a.uns.get("pca", {}).get("variance_ratio", []), dtype=float)
            n_pcs, rationale = _pick_pca_elbow(vr, max_n=n_pcs_max)
        else:
            try:
                n_pcs = int(n_pcs_param)
            except (TypeError, ValueError):
                n_pcs, rationale = _pick_pca_elbow(
                    np.asarray(a.uns.get("pca", {}).get("variance_ratio", []), dtype=float),
                    max_n=n_pcs_max)
                rationale = f"fallback to elbow ({rationale}); plan value {n_pcs_param!r} invalid"
            else:
                rationale = f"plan-fixed n_pcs={n_pcs}"

        # Trim PCA representation to the chosen n_pcs (keeps sc.pp.neighbors fast).
        if "X_pca" in a.obsm and a.obsm["X_pca"].shape[1] > n_pcs:
            a.obsm["X_pca"] = a.obsm["X_pca"][:, :n_pcs]
            if "PCs" in a.varm and a.varm["PCs"].shape[1] > n_pcs:
                a.varm["PCs"] = a.varm["PCs"][:, :n_pcs]
            if isinstance(a.uns.get("pca"), dict):
                a.uns["pca"]["variance_ratio"] = a.uns["pca"]["variance_ratio"][:n_pcs]
                a.uns["pca"]["variance"] = a.uns["pca"].get("variance",
                    np.array([]))[:n_pcs] if isinstance(a.uns["pca"].get("variance"), np.ndarray) else a.uns["pca"].get("variance")

        sc.pp.neighbors(a, n_neighbors=n_neighbors, n_pcs=n_pcs)
        _io.write_h5ad_safe(a, art / "rna_dimred.h5ad")
        _prov.set_param(params_path, "s6_dimred.rna_n_pcs", int(n_pcs),
                        source="derived", confidence="high",
                        rationale=rationale,
                        method={"name": "_pick_pca_elbow",
                                "code_ref": "executor/stages/s6_dimred.py"})
        _prov.set_param(params_path, "s6_dimred.rna_scale_applied", bool(do_scale),
                        source="derived", confidence="high",
                        rationale=("sc.pp.scale(max_value=10) applied" if do_scale
                                   else "rna_scale=False; PCA on log-normalized but unscaled data"),
                        method={"name": "scanpy.pp.scale",
                                "code_ref": "executor/stages/s6_dimred.py"})
    else:
        # atac_only — produce an empty placeholder so the Snakemake output exists.
        import scipy.sparse as sp
        _io.write_h5ad_safe(ad.AnnData(X=sp.csr_matrix((0, 0))), art / "rna_dimred.h5ad")
        n_pcs = 0

    # --- ATAC ---
    if has_atac:
        import snapatac2 as snap
        from ..atac_latent import ATAC_LATENT_KEY
        atac_h5 = run_dir / "internal" / "artifacts" / "s5_atac_spectral" / "atac_spectral.h5ad"
        if not atac_h5.exists():
            atac_h5 = run_dir / "internal" / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"
        if atac_h5.exists():
            adata = snap.read(str(atac_h5))
            try:
                # SnapATAC2 default use_rep='X_spectral' (trimmed in S5 when drop_first=True).
                snap.pp.knn(adata, n_neighbors=n_neighbors, use_rep=ATAC_LATENT_KEY)
            except Exception as e:
                log_event(run_dir, {"stage": "s6_dimred", "event": "atac_knn_failed", "error": str(e)})
            try:
                adata.close()
            except Exception:
                pass

    log_event(run_dir, {"stage": "s6_dimred", "event": "done",
                        "n_pcs": n_pcs, "n_neighbors": n_neighbors, "branch": branch})
    return {"n_pcs": n_pcs, "n_neighbors": n_neighbors, "branch": branch}
