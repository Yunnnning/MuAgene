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
    _SCRUBLET_HVG_CAP = 3000  # max genes passed to Scrublet; HVG-filter above this
    if has_rna:
        rna = ad.read_h5ad(run_dir / "internal" / "artifacts" / "s1_rna_qc" / "rna_qc.h5ad")
        raw_counts = rna.layers["counts"] if "counts" in rna.layers else rna.X
        counts = _as_sparse(raw_counts)
        # Scrublet doubles the matrix for simulated doublets: limit gene count to avoid
        # OOM on large datasets. Identify HVGs on the count matrix when n_vars > cap.
        n_hvg_used = counts.shape[1]
        if counts.shape[1] > _SCRUBLET_HVG_CAP:
            import scanpy as sc
            tmp = ad.AnnData(X=counts.copy())
            sc.pp.normalize_total(tmp, target_sum=1e4)
            sc.pp.log1p(tmp)
            sc.pp.highly_variable_genes(tmp, n_top_genes=_SCRUBLET_HVG_CAP, flavor="seurat",
                                        inplace=True)
            hvg_mask = tmp.var["highly_variable"].values
            counts = counts[:, hvg_mask]
            n_hvg_used = int(hvg_mask.sum())
            del tmp
        rate_param = plan["stages"]["s3_doublets"]["parameters"]["scrublet_expected_rate"]["value"]
        expected_rate, rate_reason = _resolve_doublet_rate(rate_param, int(rna.n_obs))
        _prov.set_param(params_path, "s3_doublets.scrublet_expected_rate_resolved",
                        float(expected_rate),
                        source="derived", confidence="high",
                        rationale=(f"Resolved from plan value={rate_param!r}: {rate_reason}. "
                                   "Tracks 10x's ~0.8%/1000 cells empirical doublet rate. "
                                   f"Scrublet run on top {n_hvg_used} HVGs."),
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

    # ---- Branch: per-modality independent (separate) vs unified policy ----
    n_dropped_rna_at_join = 0
    n_dropped_atac_at_join = 0
    n_joint: int | None = None

    if workflow_branch == "separate":
        # Independent per-modality removal.
        # RNA and ATAC are from independent samples with disjoint barcodes;
        # cross-modal union policy is undefined here.
        rna_df = pd.DataFrame({
            "barcode": rna.obs_names.to_numpy() if rna is not None else np.array([], dtype=object),
            "scrublet_score": scores,
            "scrublet_is_doublet": flags.astype(bool) if flags.size else flags,
            "removed": flags.astype(bool) if flags.size else flags,
        })
        atac_df = pd.DataFrame({
            "barcode": atac_bc,
            "atac_doublet_score": atac_scores,
            "atac_is_doublet": atac_flags.astype(bool) if atac_flags.size else atac_flags,
            "removed": atac_flags.astype(bool) if atac_flags.size else atac_flags,
        })
        n_rna_removed = int(rna_df["removed"].sum()) if has_rna else 0
        n_atac_removed = int(atac_df["removed"].sum()) if has_atac else 0
        n_removed = n_rna_removed + n_atac_removed
        chosen_policy = "independent"
        overlap = None

        # calls.parquet: concat both modalities; fill cross-modal columns with NaN.
        combined = pd.concat([
            rna_df.assign(atac_doublet_score=float("nan"), atac_is_doublet=False,
                          chosen_policy=chosen_policy),
            atac_df.assign(scrublet_score=float("nan"), scrublet_is_doublet=False,
                           chosen_policy=chosen_policy),
        ], ignore_index=True)
        combined.to_parquet(art / "calls.parquet")
        (art / "overlap_summary.json").write_text(json.dumps({
            "branch": "separate",
            "policy": "independent",
            "rationale": ("separate branch: modalities are independent samples with disjoint "
                          "barcodes; each modality's doublets are removed by its own detector."),
            "n_rna_removed": n_rna_removed,
            "n_atac_removed": n_atac_removed,
        }, indent=2))

        _prov.set_param(params_path, "s3_doublets.removal_policy", "independent",
                        source="derived", confidence="high",
                        rationale=("separate branch: modalities are from independent samples "
                                   "with disjoint barcodes. Each modality's doublets are removed "
                                   "by its own detector (Scrublet for RNA, SnapATAC2 for ATAC). "
                                   "Cross-modal union reconciliation does not apply."),
                        method={"name": "s3.independent_per_modality",
                                "code_ref": "executor/stages/s3_doublets.py"})

        rm_rna = set(rna_df.loc[rna_df["removed"], "barcode"])
        rm_atac = set(atac_df.loc[atac_df["removed"], "barcode"])
        rna_survivors = (set(rna.obs_names) - rm_rna) if (has_rna and rna is not None) else set()
        atac_survivors = (set(atac_bc) - rm_atac) if (has_atac and atac_bc) else set()

    else:
        # paired / rna_only / atac_only — unified per-barcode table + cross-modal policy.
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

        overlap = _pol.four_way_overlap(merged["scrublet_is_doublet"], merged["atac_is_doublet"])
        study_goal = plan["stages"]["s3_doublets"]["parameters"]["study_goal"]["value"]
        policy = _pol.recommend_policy(study_goal)
        chosen_policy = policy["recommendation"]
        removed = _pol.apply_policy(merged["scrublet_is_doublet"], merged["atac_is_doublet"], chosen_policy)
        merged["chosen_policy"] = chosen_policy
        merged["removed"] = removed
        n_removed = int(removed.sum())

        merged.to_parquet(art / "calls.parquet")
        (art / "overlap_summary.json").write_text(json.dumps({
            "overlap": overlap,
            "study_goal": study_goal,
            "recommended_policy": chosen_policy,
            "rationale": policy["rationale"],
            "n_removed": n_removed,
        }, indent=2))

        _prov.set_param(params_path, "s3_doublets.removal_policy", chosen_policy,
                        source="recommended", confidence="high",
                        rationale=policy["rationale"])
        _prov.set_param(params_path, "s3_doublets.overlap", overlap,
                        source="derived", confidence="high",
                        rationale=(f"Four-way overlap of Scrublet (RNA) vs {atac_method} (ATAC)"),
                        method={"name": "doublet_policy.four_way_overlap",
                                "code_ref": "executor/methods/doublet_policy.py"})

        rm_bc = set(merged.loc[merged["removed"], "barcode"])
        rna_survivors = (set(rna.obs_names) - rm_bc) if (has_rna and rna is not None) else set()
        atac_survivors = (set(atac_bc) - rm_bc) if (has_atac and atac is not None) else set()

        # Paired branch: intersect survivor sets so S4..S8 see matched barcodes.
        if workflow_branch == "paired":
            joint_bc = rna_survivors & atac_survivors
            if not joint_bc:
                raise ValueError(
                    "S3: paired-branch joint barcode intersection is empty after QC + "
                    f"doublet removal (n_rna_survivors={len(rna_survivors)}, "
                    f"n_atac_survivors={len(atac_survivors)}). Check that S0 established "
                    "real cell-level pairing (look at `ingest.pairing_decision.method` in "
                    "parameters.yaml) — if it did, your S1/S2 QC thresholds may have "
                    "filtered away the shared cells; revise via `executor revise s1_rna_qc "
                    "...` or `executor revise s2_atac_qc ...`. If S0 did NOT establish "
                    "pairing (method=pairing.translation_table required when whitelists "
                    "differ), supply `barcode_translation_path` in run.yaml and rerun S0."
                )
            n_dropped_rna_at_join = len(rna_survivors) - len(joint_bc)
            n_dropped_atac_at_join = len(atac_survivors) - len(joint_bc)
            rna_survivors = joint_bc
            atac_survivors = joint_bc
            n_joint = len(joint_bc)
            log_event(run_dir, {"stage": "s3_doublets", "event": "paired_intersection",
                                 "n_joint": n_joint,
                                 "n_dropped_rna_at_join": n_dropped_rna_at_join,
                                 "n_dropped_atac_at_join": n_dropped_atac_at_join})
            _prov.set_param(params_path, "s3_doublets.paired_intersection",
                            {"n_joint": int(n_joint),
                             "n_dropped_rna_at_join": int(n_dropped_rna_at_join),
                             "n_dropped_atac_at_join": int(n_dropped_atac_at_join)},
                            source="derived", confidence="high",
                            rationale=("Paired branch: RNA and ATAC barcodes intersected after "
                                       "doublet removal so downstream stages (S4-S8) operate on the "
                                       "joint cell set. Cells in only one modality's post-doublet "
                                       "set are dropped here rather than at S8 assembly."),
                            method={"name": "s3.paired_intersection",
                                    "code_ref": "executor/stages/s3_doublets.py"})
            (art / "joint_barcodes.txt").write_text("\n".join(sorted(joint_bc)) + ("\n" if joint_bc else ""))

    # Record ATAC detection method whenever ATAC ran (branch-independent).
    if has_atac:
        _prov.set_param(params_path, "s3_doublets.atac_method", atac_method,
                        source="derived", confidence="high",
                        rationale=("Plan named AMULET (fragment multi-allelic overlap). "
                                   "Current environment uses SnapATAC2's native scrublet + "
                                   "filter_doublets instead (Scrublet adapted to ATAC tile "
                                   "matrix). Deviation explicitly recorded; raw scores + flags "
                                   "preserved in calls.parquet."),
                        method={"name": atac_method,
                                "code_ref": "executor/stages/s3_doublets.py::_score_atac_doublets_snapatac"})

    # Apply removal to RNA (+ write filtered h5ad). Always produce the rna_post
    # sentinel (empty for atac_only) so the declared Snakemake output exists.
    rna_out = art / "rna_post_doublet.h5ad"
    n_rna_post = 0
    if has_rna and rna is not None:
        keep_rna = rna.obs_names.isin(rna_survivors)
        rna_f = rna[keep_rna].copy()
        rna_f.write_h5ad(rna_out)
        n_rna_post = int(rna_f.n_obs)
    else:
        # atac_only — write an empty placeholder AnnData so Snakemake's declared
        # output exists. Downstream stages consult workflow_branch and skip the
        # RNA path entirely.
        import scipy.sparse as sp
        ad.AnnData(X=sp.csr_matrix((0, 0))).write_h5ad(rna_out)

    # ATAC: subset via SnapATAC2 to surviving cells (only if ATAC path ran).
    atac_out = art / "atac_post_doublet.h5ad"
    if has_atac and atac is not None:
        atac_keep_idx = [i for i, bc in enumerate(atac_bc) if bc in atac_survivors]
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
        elif workflow_branch != "paired" and len(atac_survivors) == len(atac_bc):
            # No cells filtered (non-paired, no doublets removed) — copy input.
            import shutil
            atac.close()
            atac = None
            atac_in = run_dir / "internal" / "artifacts" / "s2_atac_qc" / "atac_qc.h5ad"
            shutil.copy(atac_in, atac_out)
        else:
            # Empty survivor set (extreme edge: zero joint cells). Write a tiny
            # placeholder so the declared Snakemake output exists.
            import scipy.sparse as sp
            ad.AnnData(X=sp.csr_matrix((0, 0))).write_h5ad(atac_out)
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
                         "n_removed": n_removed,
                         "n_rna_post": n_rna_post,
                         "n_joint": n_joint,
                         "branch": workflow_branch})
    return {"policy": chosen_policy, "overlap": overlap,
            "n_removed": n_removed, "n_joint": n_joint}
