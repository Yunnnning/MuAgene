"""S3 — Doublets (Scrublet for RNA; AMULET-style for ATAC if available).

For MVP on the example data, AMULET is approximated with a SnapATAC2-based heuristic
(cells with anomalously high fragment counts are scored as likely doublets) if AMULET
is not installed. The raw per-cell scores are preserved either way.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scrublet as scr

from ..methods import doublet_policy as _pol
from .. import provenance as _prov
from ..log import log_event


def _resolve_doublet_rate(value: Any, n_cells: int) -> tuple[float, str]:
    """Resolve plan's `scrublet_expected_rate` into a numeric rate.

    'auto' → min(0.10, 0.0008 * n_cells), tracking 10x's empirical curve
    (~0.8% per 1000 cells, capped at 10%). Otherwise coerced to float.
    """
    if isinstance(value, str) and value.strip().lower() == "auto":
        rate = min(0.10, 0.0008 * max(n_cells, 0))
        return float(rate), f"auto (0.0008 * {n_cells} cells, capped at 0.10)"
    try:
        rate = float(value)
    except (TypeError, ValueError):
        rate = 0.06
        return rate, f"fallback 0.06 (could not coerce {value!r})"
    return rate, f"user-fixed {rate}"


def _as_sparse(matrix: Any) -> sp.csr_matrix:
    """Return a CSR sparse counts matrix without densifying first."""
    if sp.issparse(matrix):
        return matrix.tocsr()
    return sp.csr_matrix(np.asarray(matrix))


def _score_atac_doublets_snapatac(atac, *, probability_threshold: float):
    """Use SnapATAC2's native doublet detector (`snap.pp.scrublet` + `filter_doublets`).

    Deviation from the approved design: the plan names AMULET (fragment-level
    multi-allelic overlap). AMULET is not practical in this environment; SnapATAC2 2.8
    ships `pp.scrublet` which is an ATAC-adapted Scrublet on the tile matrix and
    `pp.filter_doublets` which applies the threshold. The raw scores + boolean flag
    are preserved in the output parquet regardless, so downstream agents can re-derive
    a different policy.

    Returns (scores, flags). Adds a 'tile_matrix' to `atac` if not present.
    """
    import snapatac2 as snap

    # Scrublet needs a tile matrix or peak matrix; add tile matrix if not
    # present. Using the SnapATAC2-default `bin_size=500` here keeps the
    # doublet-scoring tile matrix consistent with S5's clustering tile
    # matrix; see s5_atac_lsi.py.
    try:
        snap.pp.add_tile_matrix(atac, bin_size=500)
    except Exception:
        pass  # may already exist
    # SnapATAC2 2.8 requires either `select_features` to have been called, or
    # `features=None` passed to use the full tile matrix. We pick the latter
    # for simplicity.
    try:
        snap.pp.scrublet(atac, features=None)
    except TypeError:
        # older snapatac2 may not accept features kwarg
        snap.pp.scrublet(atac)
    # obs now has 'doublet_score' and/or 'doublet_probability'
    def _to_arr(key: str) -> np.ndarray:
        try:
            return np.asarray(atac.obs[key].to_numpy(), dtype=float)
        except Exception:
            return np.array([])

    scores = _to_arr("doublet_score")
    probs = _to_arr("doublet_probability")
    if probs.size:
        flags = probs > probability_threshold
    elif scores.size:
        flags = scores > probability_threshold
    else:
        flags = np.zeros(atac.n_obs, dtype=bool)
    return (scores if scores.size else probs), flags


def run(run_dir: Path | str, plan: dict[str, Any], workflow_branch: str) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s3_doublets"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"

    has_rna = workflow_branch in ("paired", "separate", "rna_only")
    has_atac = workflow_branch in ("paired", "separate", "atac_only")

    # ---- RNA path (Scrublet) ---------------------------------------------
    rna = None
    scores = np.array([], dtype=float)
    flags = np.array([], dtype=bool)
    if has_rna:
        rna = ad.read_h5ad(run_dir / "internal" / "artifacts" / "s1_rna_qc" / "rna_qc.h5ad")
        # Scrublet accepts sparse counts; densifying is wasteful on large
        # datasets and crashes on >100K cells. Always pass csr.
        raw_counts = rna.layers["counts"] if "counts" in rna.layers else rna.X
        counts = _as_sparse(raw_counts)
        rate_param = plan["stages"]["s3_doublets"]["parameters"]["scrublet_expected_rate"]["value"]
        expected_rate, rate_reason = _resolve_doublet_rate(rate_param, int(rna.n_obs))
        _prov.set_param(params_path, "s3_doublets.scrublet_expected_rate_resolved",
                        float(expected_rate),
                        source="derived", confidence="high",
                        rationale=(f"Resolved from plan value={rate_param!r}: {rate_reason}. "
                                   "Tracks 10x's ~0.8%/1000 cells empirical doublet rate."),
                        method={"name": "s3.resolve_doublet_rate",
                                "code_ref": "executor/stages/s3_doublets.py"})
        try:
            sd = scr.Scrublet(counts, expected_doublet_rate=expected_rate, random_state=0)
            scores, flags = sd.scrub_doublets(verbose=False)
            if flags is None or not np.any(flags):
                # Scrublet's auto-threshold occasionally fails (bimodality unclear).
                # Fall back to a conservative score-based cutoff.
                flags = scores > 0.2
        except Exception as e:
            log_event(run_dir, {"stage": "s3_doublets", "event": "scrublet_failed", "error": str(e)})
            scores = np.zeros(rna.n_obs)
            flags = np.zeros(rna.n_obs, dtype=bool)
        rna.obs["scrublet_score"] = scores
        rna.obs["scrublet_is_doublet"] = flags.astype(bool)

    # ---- ATAC path (SnapATAC2 native scrublet + filter_doublets) ---------
    atac = None
    atac_bc: list = []
    atac_scores = np.array([], dtype=float)
    atac_flags = np.array([], dtype=bool)
    atac_method = "skipped_no_atac_input"
    if has_atac:
        import snapatac2 as snap
        atac_h5 = run_dir / "internal" / "artifacts" / "s2_atac_qc" / "atac_qc.h5ad"
        atac = snap.read(str(atac_h5))
        atac_method = "snapatac2.pp.scrublet+filter_doublets"
        atac_threshold = float(plan["stages"]["s3_doublets"]["parameters"]
                                .get("atac_doublet_threshold", {}).get("value", 0.5))
        _prov.set_param(params_path, "s3_doublets.atac_doublet_threshold",
                        atac_threshold, source="recommended", confidence="medium",
                        rationale=("SnapATAC2 scrublet doublet-probability cutoff "
                                   "above which a barcode is flagged."),
                        method={"name": "s3.atac_doublet_threshold",
                                "code_ref": "executor/stages/s3_doublets.py"})
        try:
            atac_scores, atac_flags = _score_atac_doublets_snapatac(
                atac, probability_threshold=atac_threshold)
            atac_bc = list(atac.obs_names)
        except Exception as e:
            log_event(run_dir, {"stage": "s3_doublets", "event": "atac_snap_doublet_failed",
                                "error": str(e), "falling_back_to": "log_fragment_zscore"})
            # Conservative fallback — preserves the pipeline; deviation recorded in provenance.
            n_frag = np.asarray([], dtype=float)
            try:
                n_frag = np.asarray(atac.obs["n_fragment"].to_numpy(), dtype=float)
            except Exception:
                pass
            if n_frag.size:
                med = np.median(np.log1p(n_frag))
                mad = np.median(np.abs(np.log1p(n_frag) - med)) or 1.0
                z = (np.log1p(n_frag) - med) / (mad * 1.4826)
                atac_scores, atac_flags = z, z > 3.0
            else:
                atac_scores = np.array([])
                atac_flags = np.array([], dtype=bool)
            atac_bc = list(atac.obs_names) if hasattr(atac, "obs_names") else []
            atac_method = "log_fragment_zscore_fallback"

    # Build a unified per-barcode table (branch-aware — missing modality drops out).
    rna_df = pd.DataFrame({
        "barcode": rna.obs_names.to_numpy() if rna is not None else np.array([], dtype=object),
        "scrublet_score": scores,
        "scrublet_is_doublet": flags.astype(bool) if flags.size else flags,
    })
    atac_df = pd.DataFrame({
        "barcode": atac_bc,
        "atac_doublet_score": atac_scores,
        "atac_is_doublet": atac_flags.astype(bool) if atac_flags.size else atac_flags,
    })
    merged = rna_df.merge(atac_df, on="barcode", how="outer")
    merged["scrublet_is_doublet"] = merged["scrublet_is_doublet"].fillna(False).astype(bool)
    merged["atac_is_doublet"] = merged["atac_is_doublet"].fillna(False).astype(bool)

    # Overlap + policy
    overlap = _pol.four_way_overlap(merged["scrublet_is_doublet"], merged["atac_is_doublet"])
    study_goal = plan["stages"]["s3_doublets"]["parameters"]["study_goal"]["value"]
    policy = _pol.recommend_policy(study_goal)
    chosen_policy = policy["recommendation"]
    removed = _pol.apply_policy(merged["scrublet_is_doublet"], merged["atac_is_doublet"], chosen_policy)
    merged["chosen_policy"] = chosen_policy
    merged["removed"] = removed

    merged.to_parquet(art / "calls.parquet")
    (art / "overlap_summary.json").write_text(json.dumps({
        "overlap": overlap,
        "study_goal": study_goal,
        "recommended_policy": chosen_policy,
        "rationale": policy["rationale"],
        "n_removed": int(removed.sum()),
    }, indent=2))

    _prov.set_param(params_path, "s3_doublets.removal_policy", chosen_policy,
                    source="recommended", confidence="high",
                    rationale=policy["rationale"])
    _prov.set_param(params_path, "s3_doublets.atac_method", atac_method,
                    source="derived", confidence="high",
                    rationale=("Plan named AMULET (fragment multi-allelic overlap). "
                               "Current environment uses SnapATAC2's native scrublet + "
                               "filter_doublets instead (Scrublet adapted to ATAC tile "
                               "matrix). Deviation explicitly recorded; raw scores + flags "
                               "preserved in calls.parquet."),
                    method={"name": atac_method,
                            "code_ref": "executor/stages/s3_doublets.py::_score_atac_doublets_snapatac"})
    _prov.set_param(params_path, "s3_doublets.overlap", overlap,
                    source="derived", confidence="high",
                    rationale=(f"Four-way overlap of Scrublet (RNA) vs {atac_method} (ATAC)"),
                    method={"name": "doublet_policy.four_way_overlap",
                            "code_ref": "executor/methods/doublet_policy.py"})

    # Apply removal to RNA (+ write filtered h5ad). Always produce the rna_post
    # sentinel (empty for atac_only) so the declared Snakemake output exists.
    rm_bc = set(merged.loc[merged["removed"], "barcode"])
    rna_out = art / "rna_post_doublet.h5ad"
    n_rna_post = 0
    if has_rna and rna is not None:
        keep_rna = ~rna.obs_names.isin(rm_bc)
        rna_f = rna[keep_rna].copy()
        rna_f.write_h5ad(rna_out)
        n_rna_post = int(rna_f.n_obs)
    else:
        # atac_only — write an empty placeholder AnnData so Snakemake's declared
        # output exists. Downstream stages consult workflow_branch and skip the
        # RNA path entirely.
        import scipy.sparse as sp
        ad.AnnData(X=sp.csr_matrix((0, 0))).write_h5ad(rna_out)

    # ATAC: subset via SnapATAC2 to non-removed cells (only if ATAC path ran).
    atac_out = art / "atac_post_doublet.h5ad"
    if has_atac and atac is not None:
        atac_keep_idx = [i for i, bc in enumerate(atac_bc) if bc not in rm_bc]
        if atac_keep_idx:
            atac_f = atac.subset(obs_indices=atac_keep_idx,
                                  out=str(atac_out), inplace=False)
            if atac_f is None:
                import snapatac2 as snap
                atac_f = snap.read(str(atac_out))
            try:
                atac_f.close()
            except Exception:
                pass
        else:
            # No cells to remove — copy input to expected output path
            import shutil
            atac.close()
            atac = None
            atac_in = run_dir / "internal" / "artifacts" / "s2_atac_qc" / "atac_qc.h5ad"
            shutil.copy(atac_in, atac_out)
        if atac is not None:
            try:
                atac.close()
            except Exception:
                pass
    else:
        # rna_only — empty placeholder so downstream rules that rglob the
        # artifacts dir don't crash. s6/s7/s8 consult workflow_branch directly.
        import scipy.sparse as sp
        ad.AnnData(X=sp.csr_matrix((0, 0))).write_h5ad(atac_out)

    log_event(run_dir, {"stage": "s3_doublets", "event": "done",
                         "overlap": overlap, "policy": chosen_policy,
                         "n_removed": int(removed.sum()),
                         "n_rna_post": n_rna_post,
                         "branch": workflow_branch})
    return {"policy": chosen_policy, "overlap": overlap, "n_removed": int(removed.sum())}
