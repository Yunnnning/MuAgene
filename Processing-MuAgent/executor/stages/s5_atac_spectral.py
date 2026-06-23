"""S5 — ATAC spectral embedding via SnapATAC2, plus flexible feature export.

SnapATAC2 backed AnnData writes in-place on the file it was opened from. To preserve
the S3 output unchanged and make S6 (PCA + neighbors)+ work on a distinct artifact, we COPY the S3
output to s5/atac_spectral.h5ad then modify that copy in place.

S5 performs two independent operations:

  (1) CLUSTERING LATENT: tile matrix → feature selection → snap.tl.spectral
      (Laplacian eigenmaps with IDF feature weights) → `adata.obsm['X_spectral']`.
      When `drop_first=True`, the first component is removed from X_spectral so
      SnapATAC2 defaults (knn, umap, leiden) see the trimmed embedding. A copy
      is also stored as `X_lsi` for backward compatibility.

  (2) FEATURE EXPORT: prefer a peak-by-cell matrix for downstream data
      integration. Priority order for the peak BED coordinates:
        0. User-supplied peaks (atac_peaks_path in run.yaml)
        1. Cell Ranger ARC pre-called peaks (single_file_multiome shortcut)
        2. Peaks pre-called by S2 ATAC QC for FRiP computation
           (s2_atac_qc/peaks_macs3.bed or peaks_arc.bed)
      If all peak sources are absent or fail, S5 falls back to the verified
      tile matrix that fed the spectral step. Only if that also fails does S5
      emit no feature matrix and let S8 surface a latent-only ATAC AnnData.
      S5 no longer calls MACS3 independently — peak calling is S2's responsibility.

  Outputs written to `s5_atac_spectral/`:
    feature_matrix.npz   — scipy.sparse.csr_matrix (cells × peaks or tiles).
    feature_names.tsv    — one interval per line.
    feature_kind.txt     — "peak_matrix" | "tile_matrix" | "".
"""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

from .. import io as _io
from .. import provenance as _prov
from ..atac_latent import ATAC_LATENT_ALIAS, ATAC_LATENT_KEY
from ..log import log_event


def _resolve_add_chr_prefix(params_path: Path | str) -> bool:
    """Whether S2 renamed ATAC fragments Ensembl→UCSC — peaks must match it.

    Prefers the explicit ``s2_atac_qc.add_chr_prefix`` param. For runs whose S2
    predates that param, derive it from the recorded prepared-fragments filename:
    ``io.prepare_fragments_for_snapatac`` writes the ``_chrnorm.tsv.gz`` suffix iff
    the Ensembl→UCSC rename was applied, so the suffix is an authoritative signal.
    """
    val = _prov.get_value(params_path, "s2_atac_qc.add_chr_prefix", None)
    if val is not None:
        return bool(val)
    cbf = str(_prov.get_value(params_path, "s2_atac_qc.chrom_bound_filter", "") or "")
    return cbf.endswith("_chrnorm.tsv.gz")


def _build_atac_working(run_dir: Path, dst: Path):
    """Materialise S5's working ATAC object (atac_spectral.h5ad) as a SnapATAC2-native
    backed file from the canonical post-QC h5mu's ``atac`` modality (fragments + chrom
    sizes). The caller then runs add_tile_matrix (which rebuilds the tile matrix in
    backed mode) + spectral on it, and S6/S8 reopen this same snap-native file. A
    snap-native file is required here: persisting the tile matrix via anndata is ~2x
    larger and OOMs S6's ``snap.read``.

    Legacy / transition fallback: copy the transient s3 ``atac_post_doublet.h5ad`` (or
    the s2 ``atac_qc.h5ad``) when the h5mu is absent (a pre-dedup run, or qc_handoff
    not yet run on this run dir).
    """
    import snapatac2 as snap
    from ..run_paths import RunPaths
    h5mu = RunPaths(run_dir).post_qc_h5mu
    if h5mu.exists():
        import mudata as mu
        mod = mu.read_h5ad(str(h5mu), "atac")  # fragments in obsm, chrom sizes in uns
        barcodes = list(mod.obs_names)
        sd = snap.AnnData(filename=str(dst), obs=mod.obs.copy())
        sd.obs_names = barcodes  # snap.AnnData ignores the DataFrame index; set explicitly
        for key in list(mod.obsm.keys()):  # fragment_paired / fragment_single
            sd.obsm[key] = mod.obsm[key]
        if "reference_sequences" in mod.uns:
            sd.uns["reference_sequences"] = mod.uns["reference_sequences"]
        return sd
    src = run_dir / "internal" / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"
    if not src.exists():
        src = run_dir / "internal" / "artifacts" / "s2_atac_qc" / "atac_qc.h5ad"
    shutil.copy(src, dst)
    _io.sync_path(dst)
    return snap.read(str(dst))


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    import snapatac2 as snap
    import json
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s5_atac_spectral"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"
    branch = _prov.current_branch(str(params_path))
    # Peaks must use the SAME chrom convention S2 applied to the fragments
    # (Ensembl→UCSC), else every peak silently fails to overlap.
    add_chr_prefix = _resolve_add_chr_prefix(params_path)

    if branch == "rna_only":
        _io.write_text_safe(art / "spectral_summary.json", json.dumps({
            "stage": "s5_atac_spectral",
            "skipped": True,
            "reason": "rna_only branch has no ATAC modality",
        }, indent=2))
        log_event(run_dir, {"stage": "s5_atac_spectral", "event": "skipped_no_atac",
                            "branch": branch})
        return {"skipped": True, "branch": branch}

    dst = art / "atac_spectral.h5ad"
    if dst.exists():
        dst.unlink()
    # Build the working object from the canonical post-QC h5mu (snap-native; the
    # existing add_tile_matrix call below rebuilds the tile matrix). Legacy runs with
    # the transient s3 h5ad still present fall back to the old copy path.
    adata = _build_atac_working(run_dir, dst)

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
    peak_source = ""   # "user_peaks" | "arc_h5" | "s2_peaks_macs3" | "s2_peaks_arc" | "tile_matrix_fallback"
    peak_source_h5: str | None = None

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
            # Strip comment lines + apply the fragment chrom convention (shared with S2).
            prepared_peaks = _io.prepare_peaks_for_snapatac(
                user_peaks_path, art / "_user_peaks_prepared.bed",
                add_chr_prefix=add_chr_prefix,
                log=lambda d: log_event(run_dir, {"stage": "s5_atac_spectral", **d}),
            )
            pm_out = snap.pp.make_peak_matrix(
                adata, peak_file=str(prepared_peaks), inplace=False, file=str(peak_h5),
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
    except (KeyboardInterrupt, SystemExit):
        raise
    except BaseException as e:  # incl. SnapATAC2 Rust PanicException (BaseException)
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

    # ---- Priority 2: peaks pre-called by S2 ATAC QC --------------------
    # S2 calls MACS3 (or uses user/ARC peaks) to compute FRiP and writes the
    # peak BED to its artifacts directory. Reuse those coordinates here so
    # S5 never needs to call MACS3 independently; the peak set is identical
    # regardless of whether the cell set is pre- or post-doublet-removal.
    if peak_X is None:
        s2_art = run_dir / "internal" / "artifacts" / "s2_atac_qc"
        for s2_candidate, src_label in [
            (s2_art / "peaks_macs3.bed", "s2_peaks_macs3"),
            (s2_art / "peaks_arc.bed",   "s2_peaks_arc"),
        ]:
            if not s2_candidate.exists():
                continue
            try:
                peak_h5 = art / "peak_matrix_s2peaks.h5ad"
                if peak_h5.exists():
                    peak_h5.unlink()
                # Strip comments + match fragment chrom convention (shared helper).
                # The 'not chrom.startswith("chr")' guard makes this safe for both
                # Ensembl ARC peaks and already-UCSC MACS3 peaks (no double-prefix).
                prepared_s2_peaks = _io.prepare_peaks_for_snapatac(
                    s2_candidate, art / "_s2_peaks_prepared.bed",
                    add_chr_prefix=add_chr_prefix,
                    log=lambda d: log_event(run_dir, {"stage": "s5_atac_spectral", **d}),
                )
                pm_out = snap.pp.make_peak_matrix(
                    adata, peak_file=str(prepared_s2_peaks), inplace=False, file=str(peak_h5),
                )
                peak_ad = pm_out if pm_out is not None else snap.read(str(peak_h5))

                peak_obs_names = [str(v) for v in list(peak_ad.obs_names)]
                peak_bc_set = set(peak_obs_names)
                missing = [bc for bc in s5_barcodes if bc not in peak_bc_set]
                if missing:
                    raise RuntimeError(
                        f"S2 peak matrix missing {len(missing)} of {len(s5_barcodes)} "
                        f"S5 barcodes (source={s2_candidate.name}); refusing misaligned export."
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
                        f"S2 peak matrix shape {X.shape} invalid for n_obs={n_obs_expected}."
                    )
                names = [str(v) for v in list(peak_ad.var_names)]
                if len(names) != X.shape[1]:
                    raise RuntimeError(
                        f"S2 peak var_names length {len(names)} != ncols {X.shape[1]}."
                    )

                peak_X, peak_names = X, names
                peak_source = src_label
                try:
                    peak_ad.close()
                except Exception:
                    pass
                break
            except (KeyboardInterrupt, SystemExit):
                raise
            except BaseException as e:  # incl. SnapATAC2 Rust PanicException
                log_event(run_dir, {"stage": "s5_atac_spectral",
                                     "event": "s2_peaks_reuse_skipped",
                                     "candidate": s2_candidate.name, "reason": str(e)})

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
                                 "note": "no peak source available; exporting verified tile matrix"})
        except Exception as e:
            log_event(run_dir, {"stage": "s5_atac_spectral",
                                 "event": "tile_matrix_fallback_failed",
                                 "error": str(e)})

    if peak_X is not None and peak_names is not None:
        if feature_kind is None:
            feature_kind = "peak_matrix"

        sp.save_npz(str(art / "feature_matrix.npz"), peak_X)
        _io.sync_path(art / "feature_matrix.npz")
        _io.write_text_safe(art / "feature_names.tsv", "\n".join(peak_names))
        _io.write_text_safe(art / "feature_kind.txt", feature_kind)

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
                "no peak set available (user peaks, ARC peaks, and S2 pre-called peaks "
                "all absent or failed). Downstream consumers will see "
                "uns['atac_feature_kind']='tile_matrix'. "
                f"Clustering latent {ATAC_LATENT_KEY} (tile-derived) is preserved in "
                ".obsm and is unchanged."
            )
    else:
        # All paths failed — empty sidecar → S8 emits latent_only.
        _io.write_text_safe(art / "feature_kind.txt", "")
        export_reason = (
            "Feature-level ATAC export failed on all paths "
            "(user peaks, ARC peaks, S2 pre-called peaks, tile fallback). "
            f"S8 will emit a latent_only AnnData with zero-column .X and "
            f"preserved {ATAC_LATENT_KEY}."
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
                    rationale=("user_peaks: user-supplied BED via atac_peaks_path; "
                                "arc_h5: Cell Ranger ARC pre-called peaks; "
                                "s2_peaks_macs3: MACS3 peaks pre-called by S2 for FRiP; "
                                "s2_peaks_arc: ARC-derived peaks written by S2 for FRiP; "
                                "tile_matrix_fallback: all peak sources absent/failed."),
                    method={"name": "peak_source_selector",
                            "code_ref": "executor/stages/s5_atac_spectral.py"})
    if peak_source_h5 is not None:
        _prov.set_param(params_path, "s5_atac_spectral.peak_source_h5", peak_source_h5,
                        source="derived", confidence="high",
                        rationale="10x ARC h5 from which peaks were extracted.",
                        method={"name": "io.load_atac_from_10x_h5",
                                "code_ref": "executor/io.py"})

    import json as _json
    _io.write_text_safe(art / "spectral_summary.json", _json.dumps({
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
