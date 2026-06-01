"""S5 — ATAC spectral embedding via SnapATAC2, plus flexible feature export.

SnapATAC2 backed AnnData writes in-place on the file it was opened from. To preserve
the S3 output unchanged and make S6+ work on a distinct artifact, we COPY the S3
output to s5/atac_spectral.h5ad then modify that copy in place.

S5 performs two independent operations:

  (1) CLUSTERING LATENT: tile matrix → feature selection → snap.tl.spectral
      (Laplacian eigenmaps with IDF feature weights) → `adata.obsm['X_spectral']`.
      When `drop_first=True`, the first component is removed from X_spectral so
      SnapATAC2 defaults (knn, umap, leiden) see the trimmed embedding. A copy
      is also stored as `X_lsi` for backward compatibility.

  (2) FEATURE EXPORT: prefer a peak-by-cell matrix for downstream data
      integration. The preferred path is the Cell Ranger ARC h5 (fast shortcut
      when `single_file_multiome` was detected at S0). Otherwise peaks are
      called from fragments via SnapATAC2's MACS3 integration with
      `groupby=None` (all S5 cells as one peak-calling group). If BOTH peak
      paths fail, S5 falls back to the verified tile matrix that fed the
      spectral step. Only if that fallback also fails does S5 emit no feature
      matrix and let S8 surface a latent-only ATAC AnnData.

  Outputs written to `s5_atac_spectral/`:
    feature_matrix.npz   — scipy.sparse.csr_matrix (cells × peaks or tiles).
    feature_names.tsv    — one interval per line.
    feature_kind.txt     — "peak_matrix" | "tile_matrix" | "".
"""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

from .. import provenance as _prov
from ..atac_latent import ATAC_LATENT_ALIAS, ATAC_LATENT_KEY
from ..log import log_event


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    import snapatac2 as snap
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s5_atac_spectral"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"

    src = run_dir / "internal" / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"
    if not src.exists():
        src = run_dir / "internal" / "artifacts" / "s2_atac_qc" / "atac_qc.h5ad"
    dst = art / "atac_spectral.h5ad"
    if dst.exists():
        dst.unlink()
    shutil.copy(src, dst)
    adata = snap.read(str(dst))

    p = plan["stages"]["s5_atac_spectral"]["parameters"]
    n_components = int(p["n_components"]["value"])
    max_top_peaks = int(p["max_top_peaks"]["value"])
    drop_first = bool(p["drop_first"]["value"])

    # Tile matrix (peaks-free accessibility at fixed bin size)
    try:
        snap.pp.add_tile_matrix(adata, bin_size=500)
    except Exception as e:
        log_event(run_dir, {"stage": "s5_atac_spectral", "event": "tile_matrix_failed", "error": str(e)})

    # Feature selection (note: `select_features` lives under snap.pp in 2.8; verify before use).
    sel_fn = getattr(snap.pp, "select_features", None)
    if sel_fn is not None:
        try:
            sel_fn(adata, n_features=min(max_top_peaks, adata.n_vars))
        except Exception as e:
            log_event(run_dir, {"stage": "s5_atac_spectral", "event": "select_features_failed", "error": str(e)})
    else:
        log_event(run_dir, {"stage": "s5_atac_spectral", "event": "select_features_missing",
                             "note": "snap.pp.select_features not found; using all features."})

    # Spectral embedding via SnapATAC2 (not classical TF-IDF + SVD LSI).
    snap.tl.spectral(adata, n_comps=n_components)

    try:
        import numpy as np
        if ATAC_LATENT_KEY not in adata.obsm:
            raise RuntimeError(f"snap.tl.spectral did not write obsm['{ATAC_LATENT_KEY}']")
        emb = np.asarray(adata.obsm[ATAC_LATENT_KEY])
        if drop_first and emb.shape[1] > 1:
            emb = emb[:, 1:]
        # Overwrite X_spectral so SnapATAC2 defaults (knn, umap) respect drop_first.
        adata.obsm[ATAC_LATENT_KEY] = emb
        adata.obsm[ATAC_LATENT_ALIAS] = emb
    except Exception as e:
        log_event(run_dir, {"stage": "s5_atac_spectral", "event": "spectral_postprocess_failed",
                             "error": str(e)})
        raise

    _prov.set_param(params_path, "s5_atac_spectral.n_components", n_components,
                    source="default", confidence="high",
                    rationale="Plan default; snap.tl.spectral output dimensionality")
    _prov.set_param(params_path, "s5_atac_spectral.drop_first", drop_first,
                    source="default", confidence="high",
                    rationale="First spectral component correlates with depth")

    # Flexible feature export. Prefer a peak matrix (ARC h5, then MACS3 from
    # fragments), but if peak generation fails do NOT interrupt preprocessing:
    # fall back to the verified tile matrix that fed the spectral step above.
    # Tile matrix and X_spectral remain the clustering latent either way. If even
    # the fallback cannot be exported, we leave feature_kind empty and S8 emits
    # a latent-only ATAC AnnData.
    import numpy as np
    import scipy.sparse as sp

    s5_barcodes = list(adata.obs_names)
    n_obs_expected = int(adata.n_obs)
    peak_X = None
    peak_names: list[str] | None = None
    peak_source = ""   # "user_peaks" | "arc_h5" | "macs3_from_fragments" | "tile_matrix_fallback"
    peak_source_h5: str | None = None
    macs3_failures: list[str] = []

    # ---- Priority 0: user-supplied peaks BED ----------------------------
    # When the user provides `atac_peaks_path` in run.yaml, S5 trusts it and
    # builds the peak-by-cell matrix from those intervals via SnapATAC2's
    # `make_peak_matrix`. The spectral clustering latent above is unchanged.
    try:
        import yaml as _yaml
        from ..run_paths import RunPaths as _RunPaths
        cfg = _yaml.safe_load(_RunPaths(run_dir).run_yaml.read_text()) or {}
        user_peaks = cfg.get("atac_peaks_path")
        if user_peaks:
            user_peaks_path = Path(user_peaks)
            if not user_peaks_path.exists():
                raise RuntimeError(f"atac_peaks_path={user_peaks} not found on disk.")
            peak_h5 = art / "peak_matrix_user.h5ad"
            if peak_h5.exists():
                peak_h5.unlink()
            pm_out = snap.pp.make_peak_matrix(
                adata, peak_file=str(user_peaks_path), inplace=False, file=str(peak_h5),
            )
            peak_ad = pm_out if pm_out is not None else snap.read(str(peak_h5))

            peak_obs_names = [str(v) for v in list(peak_ad.obs_names)]
            peak_bc_set = set(peak_obs_names)
            missing = [bc for bc in s5_barcodes if bc not in peak_bc_set]
            if missing:
                raise RuntimeError(
                    f"user peaks matrix missing {len(missing)} of {len(s5_barcodes)} "
                    "S5 barcodes; refusing misaligned export."
                )
            if peak_obs_names != s5_barcodes:
                bc_to_idx = {bc: i for i, bc in enumerate(peak_obs_names)}
                order = [bc_to_idx[bc] for bc in s5_barcodes]
            else:
                order = list(range(len(peak_obs_names)))
            X = peak_ad.X[:]
            if not sp.issparse(X):
                X = sp.csr_matrix(X)
            X = X.tocsr()
            if order != list(range(len(peak_obs_names))):
                X = X[order, :]
            if X.shape[0] != n_obs_expected or X.shape[1] <= 0:
                raise RuntimeError(
                    f"user peaks matrix shape {X.shape} invalid for n_obs={n_obs_expected}."
                )
            names = [str(v) for v in list(peak_ad.var_names)]
            if len(names) != X.shape[1]:
                raise RuntimeError(
                    f"user peaks var_names length {len(names)} != ncols {X.shape[1]}."
                )
            peak_X, peak_names = X, names
            peak_source = "user_peaks"
            peak_source_h5 = str(user_peaks_path)
            try:
                peak_ad.close()
            except Exception:
                pass
    except Exception as e:
        log_event(run_dir, {"stage": "s5_atac_spectral",
                             "event": "user_peaks_path_skipped", "reason": str(e)})

    # ---- Priority 1: Cell Ranger ARC h5 peak matrix ---------------------
    if peak_X is None:
        try:
            import json as _json
            s0_report = _json.loads(
                (run_dir / "internal" / "artifacts" / "s0_ingest" / "validation_report.json").read_text()
            )
            if s0_report.get("single_file_multiome"):
                import yaml as _yaml
                from ..run_paths import RunPaths as _RunPaths
                cfg = _yaml.safe_load(_RunPaths(run_dir).run_yaml.read_text()) or {}
                rna_path = cfg.get("rna_path")
                if rna_path:
                    from .. import io as _io
                    peak_ad = _io.load_atac_from_10x_h5(rna_path)
                    peak_bc_set = set(peak_ad.obs_names)
                    missing = [bc for bc in s5_barcodes if bc not in peak_bc_set]
                    if missing:
                        raise RuntimeError(
                            f"{len(missing)} of {len(s5_barcodes)} S5 barcodes not found "
                            "in 10x ARC peak matrix; refusing misaligned export."
                        )
                    peak_ad = peak_ad[s5_barcodes].copy()
                    X = peak_ad.X
                    if not sp.issparse(X):
                        X = sp.csr_matrix(X)
                    X = X.tocsr()
                    if X.shape[0] != n_obs_expected or X.shape[1] <= 0:
                        raise RuntimeError(
                            f"peak matrix shape {X.shape} invalid for n_obs={n_obs_expected}."
                        )
                    names = [str(v) for v in list(peak_ad.var_names)]
                    if len(names) != X.shape[1]:
                        raise RuntimeError(
                            f"peak var_names length {len(names)} != ncols {X.shape[1]}."
                        )
                    peak_X, peak_names = X, names
                    peak_source, peak_source_h5 = "arc_h5", str(rna_path)
        except Exception as e:
            log_event(run_dir, {"stage": "s5_atac_spectral",
                                 "event": "arc_peak_path_skipped", "reason": str(e)})

    # ---- Priority 2: MACS3 peak calling from fragments ------------------
    if peak_X is None:
        try:
            genome = _prov.get_value(params_path, "ingest.genome_assembly", None)
            if not genome:
                raise RuntimeError("ingest.genome_assembly not set; cannot call MACS3 peaks.")
            genome_ref = getattr(snap.genome, genome, None)
            if genome_ref is None:
                raise RuntimeError(f"genome {genome!r} not supported by SnapATAC2.")

            macs_tempdir = art / "macs3_tmp"
            macs_tempdir.mkdir(parents=True, exist_ok=True)

            # Call peaks on all cells as one group (no clustering required at S5).
            # SnapATAC2 2.8: with groupby=None it returns a single polars DataFrame;
            # with a groupby column it returns dict[str, DataFrame]. Handle both.
            peaks_out = snap.tl.macs3(
                adata, groupby=None, inplace=False,
                qvalue=0.05, shift=-100, extsize=200,
                tempdir=macs_tempdir, n_jobs=1,
            )
            if peaks_out is None:
                raise RuntimeError("snap.tl.macs3 returned no peaks.")

            if isinstance(peaks_out, dict):
                if not peaks_out:
                    raise RuntimeError("snap.tl.macs3 returned an empty peak-set dict.")
                merged = snap.tl.merge_peaks(peaks_out, chrom_sizes=genome_ref)
            else:
                # Single-group result — already a polars DataFrame of peaks.
                merged = peaks_out

            # Write a BED file for make_peak_matrix. SnapATAC2's peak DataFrames
            # use lowercase `chrom/start/end` (MACS3-style) in some paths and
            # TitleCase in others — pick whichever the schema actually has.
            cols = [c.lower() for c in merged.columns]
            def _pick(name: str) -> str:
                i = cols.index(name)
                return merged.columns[i]
            bed_df = merged.select([_pick("chrom"), _pick("start"), _pick("end")])
            if bed_df.height == 0:
                raise RuntimeError("MACS3 produced zero peaks; refusing export.")

            peak_bed = art / "peaks_macs3.bed"
            bed_df.write_csv(str(peak_bed), separator="\t", include_header=False)

            # Build the peak × cell matrix directly from fragments.
            peak_h5 = art / "peak_matrix_snap.h5ad"
            if peak_h5.exists():
                peak_h5.unlink()
            pm_out = snap.pp.make_peak_matrix(
                adata, peak_file=str(peak_bed), inplace=False, file=str(peak_h5),
            )
            peak_ad = pm_out if pm_out is not None else snap.read(str(peak_h5))

            peak_obs_names = [str(v) for v in list(peak_ad.obs_names)]
            peak_bc_set = set(peak_obs_names)
            missing = [bc for bc in s5_barcodes if bc not in peak_bc_set]
            if missing:
                raise RuntimeError(
                    f"MACS3 peak matrix missing {len(missing)} of {len(s5_barcodes)} "
                    "S5 barcodes; refusing misaligned export."
                )

            if peak_obs_names != s5_barcodes:
                bc_to_idx = {bc: i for i, bc in enumerate(peak_obs_names)}
                order = [bc_to_idx[bc] for bc in s5_barcodes]
            else:
                order = list(range(len(peak_obs_names)))

            X = peak_ad.X[:]
            if not sp.issparse(X):
                X = sp.csr_matrix(X)
            X = X.tocsr()
            if order != list(range(len(peak_obs_names))):
                X = X[order, :]
            aligned_barcodes = [peak_obs_names[i] for i in order]
            if aligned_barcodes != s5_barcodes:
                raise RuntimeError(
                    "MACS3 peak matrix did not align to S5 barcodes after reorder."
                )
            if X.shape[0] != n_obs_expected or X.shape[1] <= 0:
                raise RuntimeError(
                    f"MACS3 peak matrix shape {X.shape} invalid for n_obs={n_obs_expected}."
                )
            names = [str(v) for v in list(peak_ad.var_names)]
            if len(names) != X.shape[1]:
                raise RuntimeError(
                    f"MACS3 peak var_names length {len(names)} != ncols {X.shape[1]}."
                )

            peak_X, peak_names = X, names
            peak_source = "macs3_from_fragments"
            try:
                peak_ad.close()
            except Exception:
                pass
        except Exception as e:
            macs3_failures.append(str(e))
            log_event(run_dir, {"stage": "s5_atac_spectral",
                                 "event": "macs3_peak_path_failed", "reason": str(e)})

    # ---- Sanity: spectral must have computed X_spectral ----------------
    if ATAC_LATENT_KEY not in adata.obsm:
        raise RuntimeError(
            f"{ATAC_LATENT_KEY} missing from adata.obsm — the tile/spectral clustering "
            "path failed. Refusing to emit an ATAC output without a clustering latent."
        )

    feature_kind: str | None = None
    feature_exported = False
    n_features = 0
    export_reason = ""

    # ---- Priority 3: tile-matrix fallback ------------------------------
    if peak_X is None or peak_names is None:
        try:
            n_vars_expected = int(adata.n_vars)
            if n_vars_expected <= 0:
                raise RuntimeError(
                    f"adata.n_vars={n_vars_expected} — no tile features available for fallback."
                )
            X = adata.X[:]
            if not sp.issparse(X):
                X = sp.csr_matrix(X)
            X = X.tocsr()
            if X.shape != (n_obs_expected, n_vars_expected):
                raise RuntimeError(
                    f"tile matrix shape {X.shape} does not match "
                    f"(n_obs={n_obs_expected}, n_vars={n_vars_expected})."
                )
            names = [str(v) for v in list(adata.var_names)]
            if len(names) != n_vars_expected:
                raise RuntimeError(
                    f"tile var_names length {len(names)} != ncols {n_vars_expected}."
                )

            peak_X, peak_names = X, names
            peak_source = "tile_matrix_fallback"
            feature_kind = "tile_matrix"
            log_event(run_dir, {"stage": "s5_atac_spectral",
                                 "event": "tile_matrix_fallback_engaged",
                                 "note": "peak generation failed; exporting verified tile matrix",
                                 "macs3_errors": macs3_failures})
        except Exception as e:
            log_event(run_dir, {"stage": "s5_atac_spectral",
                                 "event": "tile_matrix_fallback_failed",
                                 "error": str(e), "macs3_errors": macs3_failures})

    if peak_X is not None and peak_names is not None:
        if feature_kind is None:
            feature_kind = "peak_matrix"

        sp.save_npz(str(art / "feature_matrix.npz"), peak_X)
        (art / "feature_names.tsv").write_text("\n".join(peak_names))
        (art / "feature_kind.txt").write_text(feature_kind)

        n_features = int(peak_X.shape[1])
        feature_exported = True
        if feature_kind == "peak_matrix":
            export_reason = (
                f"Exported peak_matrix ({n_obs_expected}×{n_features}, source={peak_source}"
                f"{f', h5={peak_source_h5}' if peak_source_h5 else ''}) to "
                "s5_atac_spectral/{feature_matrix.npz, feature_names.tsv}. "
                f"Clustering latent {ATAC_LATENT_KEY} (tile-derived) is preserved in "
                ".obsm and is unchanged."
            )
        else:
            export_reason = (
                f"Exported tile_matrix fallback ({n_obs_expected}×{n_features}) — "
                f"peak generation failed (macs3_errors={macs3_failures}). Downstream "
                "consumers will see uns['atac_feature_kind']='tile_matrix'. "
                f"Clustering latent {ATAC_LATENT_KEY} (tile-derived) is preserved in "
                ".obsm and is unchanged."
            )
    else:
        # All three paths failed — empty sidecar → S8 emits latent_only.
        (art / "feature_kind.txt").write_text("")
        export_reason = (
            "Feature-level ATAC export failed on all paths (ARC, MACS3, tile fallback). "
            f"macs3_errors={macs3_failures}. S8 will emit a latent_only AnnData with "
            f"zero-column .X and preserved {ATAC_LATENT_KEY}."
        )
        feature_kind = ""

    _prov.set_param(params_path, "s5_atac_spectral.feature_matrix_exported", feature_exported,
                    source="derived", confidence="high",
                    rationale=export_reason,
                    method={"name": f"s5.feature_export.{peak_source or 'none'}",
                            "code_ref": "executor/stages/s5_atac_spectral.py"})
    _prov.set_param(params_path, "s5_atac_spectral.feature_kind", feature_kind or "",
                    source="derived", confidence="high",
                    rationale=("peak_matrix | tile_matrix | '' (latent-only) — the actual "
                               "feature representation exported by S5."),
                    method={"name": "literal", "code_ref": "executor/stages/s5_atac_spectral.py"})
    _prov.set_param(params_path, "s5_atac_spectral.peak_source", peak_source,
                    source="derived", confidence="high",
                    rationale=("arc_h5: Cell Ranger ARC pre-called peaks; "
                                "macs3_from_fragments: MACS3 called globally on S5 cells; "
                                "tile_matrix_fallback: both peak paths failed, fell back to "
                                "the verified tile matrix."),
                    method={"name": "peak_source_selector",
                            "code_ref": "executor/stages/s5_atac_spectral.py"})
    if peak_source_h5 is not None:
        _prov.set_param(params_path, "s5_atac_spectral.peak_source_h5", peak_source_h5,
                        source="derived", confidence="high",
                        rationale="10x ARC h5 from which peaks were extracted.",
                        method={"name": "io.load_atac_from_10x_h5",
                                "code_ref": "executor/io.py"})

    import json as _json
    (art / "spectral_summary.json").write_text(_json.dumps({
        "method": "snap.tl.spectral",
        "embedding_key": ATAC_LATENT_KEY,
        "embedding_alias": ATAC_LATENT_ALIAS,
        "n_components": int(n_components),
        "drop_first": bool(drop_first),
        "n_cells": int(adata.n_obs),
        "feature_matrix_exported": bool(feature_exported),
        "feature_kind": feature_kind,
        "n_features_selected": int(n_features),
        "peak_source": peak_source,
        "peak_source_h5": peak_source_h5,
    }, indent=2))
    try:
        adata.close()
    except Exception:
        pass
    log_event(run_dir, {"stage": "s5_atac_spectral", "event": "done",
                         "n_components": n_components, "drop_first": drop_first,
                         "feature_matrix_exported": feature_exported,
                         "n_features_selected": n_features})
    return {"n_components": n_components, "drop_first": drop_first,
            "feature_matrix_exported": feature_exported}
