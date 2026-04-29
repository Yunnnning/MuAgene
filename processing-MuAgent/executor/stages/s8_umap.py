"""S8 — UMAP per modality + final h5mu (paired) or separate h5ad (separate). Hard stop."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import scanpy as sc

from .. import provenance as _prov
from ..log import log_event


def load_exported_atac_features(
    feat_path: Path | str,
    names_path: Path | str,
    *,
    n_obs_expected: int,
) -> tuple[Any, list[str]] | tuple[None, None]:
    """Load S5's exported ATAC feature matrix + names, validating alignment.

    Returns `(sparse_csr_matrix, list_of_var_names)` only when:
      - both files exist,
      - the .npz loads cleanly as a sparse matrix,
      - rows equal `n_obs_expected` (row alignment),
      - columns are > 0 (no degenerate empty export),
      - names count EXACTLY equals column count (no silent trim/pad).

    Any failure returns `(None, None)` — callers must treat this as
    "no feature-level ATAC representation available" and produce an honest
    empty .X (zero columns), not a fabricated 1-column placeholder.

    Duplicate names are de-duplicated with a `-N` suffix so the resulting
    `var_names` Index is unique. MuData refuses to build a clean string-typed
    `varmap` when either modality has duplicate var_names, and `mdata.write()`
    then crashes inside `write_vlen_string_array` while serialising `varmap`.
    Deduping here keeps the fix local to the S8 consumer of S5's exports.
    """
    import scipy.sparse as sp

    feat_path = Path(feat_path)
    names_path = Path(names_path)
    if not (feat_path.exists() and names_path.exists()):
        return None, None
    try:
        mat = sp.load_npz(str(feat_path)).tocsr()
    except Exception:
        return None, None
    if mat.shape[0] != n_obs_expected:
        return None, None
    if mat.shape[1] <= 0:
        return None, None
    names = names_path.read_text().splitlines()
    if len(names) != mat.shape[1]:
        return None, None

    seen: dict[str, int] = {}
    deduped: list[str] = []
    for n in names:
        c = seen.get(n, 0)
        deduped.append(n if c == 0 else f"{n}-{c}")
        seen[n] = c + 1
    return mat, deduped


def load_atac_feature_kind(kind_path: Path | str, *, default: str = "unknown_feature_matrix") -> str:
    """Read the `feature_kind.txt` sidecar written by S5.

    S5 writes "peak_matrix" or "tile_matrix" on successful feature-matrix
    export and writes an empty sidecar when no export succeeded. `default`
    is only a safety fallback for missing / empty sidecars; it should remain
    an explicit "unknown" marker so callers never silently relabel a real load
    as peak- or tile-based. Callers must pair this with
    `load_exported_atac_features` so the kind is only surfaced to `AnnData.uns`
    when a matrix actually loaded. The latent-only branch uses the literal
    "latent_only" at the call site.
    Returns `default` when the file is missing or empty.
    """
    p = Path(kind_path)
    if not p.exists():
        return default
    s = p.read_text().strip()
    return s if s else default


def run(run_dir: Path | str, plan: dict[str, Any], workflow_branch: str) -> dict[str, Any]:
    from ..run_paths import RunPaths
    run_dir = Path(run_dir)
    paths = RunPaths(run_dir)
    art = paths.stage_dir("s8_umap")       # internal/artifacts/s8_umap/ — sentinel lives here
    art.mkdir(parents=True, exist_ok=True)
    figures = paths.deliv_figures
    figures.mkdir(parents=True, exist_ok=True)
    processed_dir = paths.deliv_processed
    processed_dir.mkdir(parents=True, exist_ok=True)
    params_path = paths.parameters_yaml

    p = plan["stages"]["s8_umap"]["parameters"]
    min_dist = float(p["min_dist"]["value"])
    spread = float(p["spread"]["value"])
    seed = int(p["random_state"]["value"])

    has_rna = workflow_branch in ("paired", "separate", "rna_only")
    has_atac = workflow_branch in ("paired", "separate", "atac_only")

    from .. import figures as _fig

    # --- RNA ---
    rna = None
    if has_rna:
        rna = ad.read_h5ad(run_dir / "internal" / "artifacts" / "s7_clustering" / "rna_clustered.h5ad")
        sc.tl.umap(rna, min_dist=min_dist, spread=spread, random_state=seed)
        if "X_umap" in rna.obsm:
            rna.obsm["X_umap_rna"] = rna.obsm["X_umap"]
        try:
            _fig.plot_umap(rna.obsm["X_umap_rna"], rna.obs["leiden_rna"],
                           title="RNA UMAP — Leiden clusters",
                           out_dir=figures, stem="s8_umap_rna_by_leiden",
                           label_name="leiden_rna")
        except Exception as e:
            log_event(run_dir, {"stage": "s8_umap", "event": "rna_plot_failed", "error": str(e)})

    # --- ATAC ---
    atac_adata = None
    if has_atac:
        import snapatac2 as snap
        atac_h5 = run_dir / "internal" / "artifacts" / "s5_atac_lsi" / "atac_lsi.h5ad"
        if not atac_h5.exists():
            atac_h5 = run_dir / "internal" / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"
        atac = snap.read(str(atac_h5))
        try:
            snap.tl.umap(atac, min_dist=min_dist, random_state=seed)
        except Exception as e:
            log_event(run_dir, {"stage": "s8_umap", "event": "atac_umap_failed", "error": str(e)})

        atac_obs_df = atac.obs[:].to_pandas()
        try:
            if "X_umap" in atac.obsm and "leiden_atac" in list(atac_obs_df.columns):
                _fig.plot_umap(np.asarray(atac.obsm["X_umap"]),
                               atac_obs_df["leiden_atac"].astype(str).to_numpy(),
                               title="ATAC UMAP — Leiden clusters",
                               out_dir=figures, stem="s8_umap_atac_by_leiden",
                               label_name="leiden_atac")
        except Exception as e:
            log_event(run_dir, {"stage": "s8_umap", "event": "atac_plot_failed", "error": str(e)})
        try:
            atac_umap = np.asarray(atac.obsm["X_umap"]) if "X_umap" in atac.obsm else None
        except Exception:
            atac_umap = None
        try:
            atac_lsi = np.asarray(atac.obsm["X_lsi"]) if "X_lsi" in atac.obsm else None
        except Exception:
            atac_lsi = None
        atac_barcodes = list(atac.obs_names)
        try:
            atac.close()
        except Exception:
            pass

        import scipy.sparse as sp
        atac_obs_df.index = atac_barcodes

        # Feature-level ATAC representation exported by S5. `load_exported_atac_features`
        # returns (None, None) if anything is wrong (missing / shape mismatch / name-count
        # mismatch / bad npz). On (None, None) we build an honest latent-only AnnData
        # with an explicit zero-column .X — we do NOT fabricate a 1-column placeholder
        # that could be mistaken for a real feature. `.obsm['X_lsi']` carries the
        # spectral latent either way.
        s5_dir = run_dir / "internal" / "artifacts" / "s5_atac_lsi"
        feat_path = s5_dir / "feature_matrix.npz"
        names_path = s5_dir / "feature_names.tsv"
        kind_path = s5_dir / "feature_kind.txt"
        atac_X, atac_var_names = load_exported_atac_features(
            feat_path, names_path, n_obs_expected=len(atac_barcodes),
        )
        latent_only = atac_X is None
        if latent_only:
            log_event(run_dir, {"stage": "s8_umap",
                                 "event": "atac_feature_matrix_unavailable",
                                 "note": "emitting latent-only AnnData with zero-column .X"})
            atac_X = sp.csr_matrix((len(atac_barcodes), 0))

        atac_adata = ad.AnnData(X=atac_X, obs=atac_obs_df)
        atac_adata.obs_names = atac_barcodes
        if not latent_only and atac_var_names is not None:
            atac_adata.var_names = atac_var_names
            # Read the kind S5 actually wrote (peak_matrix or tile_matrix).
            atac_adata.uns["atac_feature_kind"] = load_atac_feature_kind(kind_path)
        else:
            # Latent-only: mark it plainly so downstream code can branch on this.
            atac_adata.uns["atac_feature_kind"] = "latent_only"
        if atac_umap is not None:
            atac_adata.obsm["X_umap_atac"] = atac_umap
        if atac_lsi is not None:
            atac_adata.obsm["X_lsi"] = atac_lsi

    # --- Write branch-specific final output -----------------------------
    if workflow_branch == "paired":
        common = sorted(set(rna.obs_names) & set(atac_adata.obs_names))
        if common:
            rna = rna[common].copy()
            atac_adata = atac_adata[common].copy()
        import mudata as mu
        mdata = mu.MuData({"rna": rna, "atac": atac_adata})
        mdata.write(str(paths.processed_h5mu))
        outpath = str(paths.processed_h5mu)
    elif workflow_branch == "separate":
        rna.write_h5ad(paths.rna_processed_h5ad)
        atac_adata.write_h5ad(paths.atac_processed_h5ad)
        outpath = f"{paths.rna_processed_h5ad},{paths.atac_processed_h5ad}"
    elif workflow_branch == "rna_only":
        rna.write_h5ad(paths.rna_processed_h5ad)
        outpath = str(paths.rna_processed_h5ad)
    elif workflow_branch == "atac_only":
        atac_adata.write_h5ad(paths.atac_processed_h5ad)
        outpath = str(paths.atac_processed_h5ad)
    else:
        raise ValueError(f"S8: unknown workflow_branch={workflow_branch!r}")

    _prov.set_param(params_path, "s8_umap.random_state", seed, source="user",
                    confidence="high", rationale="Run seed from config")

    log_event(run_dir, {"stage": "s8_umap", "event": "done", "output": outpath,
                         "branch": workflow_branch})
    return {"output": outpath, "branch": workflow_branch}
