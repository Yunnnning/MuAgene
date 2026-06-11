"""S8 — UMAP per modality + final h5mu (paired) or separate h5ad (separate). Hard stop."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from .. import io as _io
from .. import provenance as _prov
from ..atac_latent import ATAC_LATENT_ALIAS, ATAC_LATENT_KEY, get_atac_latent
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
    processed_dir = paths.deliv_results
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
        atac_h5 = run_dir / "internal" / "artifacts" / "s5_atac_spectral" / "atac_spectral.h5ad"
        if not atac_h5.exists():
            atac_h5 = run_dir / "internal" / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"
        atac = snap.read(str(atac_h5))
        try:
            snap.tl.umap(atac, min_dist=min_dist, random_state=seed, use_rep=ATAC_LATENT_KEY)
        except Exception as e:
            log_event(run_dir, {"stage": "s8_umap", "event": "atac_umap_failed", "error": str(e)})

        atac_obs_df = atac.obs[:].to_pandas()
        atac_barcodes = list(atac.obs_names)
        atac_obs_df.index = atac_barcodes
        leiden_sidecar = (run_dir / "internal" / "artifacts" / "s7_clustering"
                          / "atac_leiden_labels.parquet")
        if "leiden_atac" not in atac_obs_df.columns and leiden_sidecar.exists():
            leiden_df = pd.read_parquet(leiden_sidecar).set_index("barcode")
            atac_obs_df = atac_obs_df.join(leiden_df, how="left")
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
            atac_emb = get_atac_latent(atac.obsm)
        except Exception:
            atac_emb = None
        try:
            atac.close()
        except Exception:
            pass

        import scipy.sparse as sp

        # Feature-level ATAC representation exported by S5. `load_exported_atac_features`
        # returns (None, None) if anything is wrong (missing / shape mismatch / name-count
        # mismatch / bad npz). On (None, None) we build an honest latent-only AnnData
        # with an explicit zero-column .X — we do NOT fabricate a 1-column placeholder
        # that could be mistaken for a real feature. `.obsm['X_spectral']` carries the
        # spectral latent either way (with `X_lsi` as a backward-compat alias).
        s5_dir = run_dir / "internal" / "artifacts" / "s5_atac_spectral"
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
        if atac_emb is not None:
            atac_adata.obsm[ATAC_LATENT_KEY] = atac_emb
            atac_adata.obsm[ATAC_LATENT_ALIAS] = atac_emb

    # --- Optional generic per-cell metadata join (any branch) -----------
    # User-supplied `cell_metadata_path` (TSV with a `barcode` column + arbitrary
    # extra columns) is left-joined into the modality .obs frames just before the
    # final write. Cells without a metadata row keep NaN. This is intentionally
    # generic — it is NOT HTO demultiplexing; the upstream pipeline must already
    # have demultiplexed if a sample_id column is expected.
    import yaml as _yaml
    cfg = _yaml.safe_load(paths.run_yaml.read_text()) or {}
    cell_metadata_path = cfg.get("cell_metadata_path")
    if cell_metadata_path:
        mpath = Path(cell_metadata_path)
        if mpath.exists():
            try:
                import pandas as _pd
                meta_df = _pd.read_csv(mpath, sep="\t", dtype=str, keep_default_na=False)
                cols_l = {c.lower(): c for c in meta_df.columns}
                bc_col = cols_l.get("barcode") or meta_df.columns[0]
                meta_df = meta_df.set_index(bc_col)
                if rna is not None:
                    rna.obs = rna.obs.join(meta_df, how="left")
                if atac_adata is not None:
                    atac_adata.obs = atac_adata.obs.join(meta_df, how="left")
                log_event(run_dir, {"stage": "s8_umap",
                                     "event": "cell_metadata_joined",
                                     "path": str(mpath),
                                     "n_columns_added": int(meta_df.shape[1])})
            except Exception as e:
                log_event(run_dir, {"stage": "s8_umap",
                                     "event": "cell_metadata_join_failed",
                                     "path": str(mpath), "error": str(e)})
        else:
            log_event(run_dir, {"stage": "s8_umap",
                                 "event": "cell_metadata_path_missing",
                                 "path": str(mpath)})

    # --- Write branch-specific final output -----------------------------
    if workflow_branch == "paired":
        # S3 intersects RNA/ATAC barcodes on the paired branch, so by here the
        # two modalities should already be aligned. The defensive check below
        # runs BEFORE MuData assembly so a divergence never reaches mdata.write,
        # and an empty intersection is a hard error with explicit context.
        rna_bcs = set(rna.obs_names)
        atac_bcs = set(atac_adata.obs_names)
        common = rna_bcs & atac_bcs
        if not common:
            raise ValueError(
                f"S8: paired output requires shared barcodes between RNA "
                f"(n={len(rna_bcs)}) and ATAC (n={len(atac_bcs)}), but the "
                "intersection is empty. This typically means S3's joint-barcode "
                "intersection was skipped (committed branch != 'paired') or that "
                "subsequent stages (S4..S7) dropped cells asymmetrically. Inspect "
                "`ingest.pairing_decision` in parameters.yaml and the s3_doublets "
                "calls.parquet / joint_barcodes.txt to diagnose."
            )
        if rna_bcs != atac_bcs:
            common_sorted = sorted(common)
            log_event(run_dir, {"stage": "s8_umap", "event": "barcode_realignment",
                                 "note": ("RNA/ATAC barcodes diverged between S3 intersection "
                                          "and S8 — falling back to intersection at assembly."),
                                 "n_rna_pre": len(rna_bcs),
                                 "n_atac_pre": len(atac_bcs),
                                 "n_common": len(common_sorted)})
            rna = rna[common_sorted].copy()
            atac_adata = atac_adata[common_sorted].copy()
        import mudata as mu
        mdata = mu.MuData({"rna": rna, "atac": atac_adata})
        _io.write_mudata_safe(mdata, paths.processed_h5mu)
        outpath = str(paths.processed_h5mu)
    elif workflow_branch == "separate":
        _io.write_h5ad_safe(rna, paths.rna_processed_h5ad)
        _io.write_h5ad_safe(atac_adata, paths.atac_processed_h5ad)
        outpath = f"{paths.rna_processed_h5ad},{paths.atac_processed_h5ad}"
    elif workflow_branch == "rna_only":
        _io.write_h5ad_safe(rna, paths.rna_processed_h5ad)
        outpath = str(paths.rna_processed_h5ad)
    elif workflow_branch == "atac_only":
        _io.write_h5ad_safe(atac_adata, paths.atac_processed_h5ad)
        outpath = str(paths.atac_processed_h5ad)
    else:
        raise ValueError(f"S8: unknown workflow_branch={workflow_branch!r}")

    _prov.set_param(params_path, "s8_umap.random_state", seed, source="user",
                    confidence="high", rationale="Run seed from config")

    log_event(run_dir, {"stage": "s8_umap", "event": "done", "output": outpath,
                         "branch": workflow_branch})
    return {"output": outpath, "branch": workflow_branch}
