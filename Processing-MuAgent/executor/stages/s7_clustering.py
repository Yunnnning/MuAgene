"""S7 — Leiden clustering with resolution sweep + stable-region knee (per modality).

This stage is split into two entry points:
  - `propose(run_dir, plan)`: runs the full resolution sweep on RNA and ATAC,
    writes sweep.parquet + sweep figures + a recommendation YAML. Writes the
    recommendation into parameters.yaml under keys `s7_clustering.rna.resolution`
    and `s7_clustering.atac.resolution` with source=derived. This does NOT
    assign final cluster labels.
  - `execute(run_dir, plan)`: reads the approved resolutions from parameters.yaml
    (which may have been revised by the user via `revise`), runs Leiden ONCE per
    modality at the approved resolution, and writes rna_clustered.h5ad with
    final `leiden_rna` / `leiden_atac` labels.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from ..methods import resolution_sweep as _rs
from .. import io as _io
from .. import provenance as _prov
from .. import resolution_compare as _rcmp
from ..log import log_event


# ---------------------------------------------------------------------------
# Propose: full sweep + recommendation, no final labels
# ---------------------------------------------------------------------------

def propose(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s7_clustering"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"
    branch = _prov.current_branch(str(params_path))
    has_rna = branch in ("paired", "separate", "rna_only")
    has_atac = branch in ("paired", "separate", "atac_only")

    p = plan["stages"]["s7_clustering"]["parameters"]
    seeds = list(p["seeds"]["value"])
    stability_floor = float(p["stability_floor"]["value"])
    rna_tilt = p["rna_tilt"]["value"]
    atac_tilt = p["atac_tilt"]["value"]
    # Per-modality grids with backwards-compatible fallback to the legacy unified grid.
    legacy_grid = list(p.get("leiden_resolution_grid", {}).get("value", [0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2]))
    rna_grid = list(p.get("leiden_resolution_grid_rna", {}).get("value", legacy_grid))
    atac_grid = list(p.get("leiden_resolution_grid_atac", {}).get("value", legacy_grid))

    # --- RNA sweep ---
    rna_rows: list[dict[str, Any]] = []
    rna_pick = {"resolution": None, "rationale": "No RNA modality in this branch", "stable_region": []}
    if has_rna:
        rna = ad.read_h5ad(run_dir / "internal" / "artifacts" / "s6_dimred" / "rna_dimred.h5ad")
        rna_rows = _rs.sweep_leiden(rna, grid=rna_grid, seeds=seeds, latent_key="X_pca")
        rna_pick = _rs.pick_stable_knee(rna_rows, tilt=rna_tilt, stability_floor=stability_floor)

    # --- ATAC sweep ---
    atac_rows: list[dict[str, Any]] = []
    atac_pick = {"resolution": None, "rationale": "No ATAC modality in this branch", "stable_region": []}
    if has_atac:
        try:
            import snapatac2 as snap
            from sklearn.metrics import adjusted_rand_score, silhouette_score
            from ..atac_latent import get_atac_latent
            atac_h5 = run_dir / "internal" / "artifacts" / "s5_atac_spectral" / "atac_spectral.h5ad"
            if not atac_h5.exists():
                atac_h5 = run_dir / "internal" / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"
            adata = snap.read(str(atac_h5))
            # Pull the spectral latent once for silhouette; subsample consistent with RNA path.
            atac_latent = None
            try:
                atac_latent = get_atac_latent(adata.obsm)
            except Exception:
                atac_latent = None
            rng = np.random.default_rng(0)
            atac_n = adata.n_obs
            sil_idx = rng.choice(atac_n, size=min(atac_n, 5000), replace=False) if atac_latent is not None else None

            assignments_per_res: dict[float, list] = {}
            for res in atac_grid:
                assignments_per_res[res] = []
                for seed in seeds:
                    snap.tl.leiden(adata, resolution=res, random_state=seed,
                                    key_added=f"_tmp_leiden_{seed}")
                    try:
                        col = adata.obs[f"_tmp_leiden_{seed}"].to_numpy()
                    except Exception:
                        col = np.asarray(adata.obs[f"_tmp_leiden_{seed}"])
                    assignments_per_res[res].append([str(x) for x in col])
            for res in atac_grid:
                base = assignments_per_res[res][0]
                aris = []
                for i in range(len(seeds)):
                    for j in range(i + 1, len(seeds)):
                        aris.append(adjusted_rand_score(assignments_per_res[res][i],
                                                        assignments_per_res[res][j]))
                # Silhouette on spectral latent (same semantics as RNA's silhouette on X_pca)
                sil = float("nan")
                if atac_latent is not None and sil_idx is not None and len(set(base)) > 1:
                    try:
                        base_arr = np.asarray(base)
                        sil = float(silhouette_score(atac_latent[sil_idx], base_arr[sil_idx]))
                    except Exception:
                        sil = float("nan")
                atac_rows.append({
                    "resolution": float(res),
                    "n_clusters": int(len(set(base))),
                    "silhouette": sil,
                    "seed_stability_ari": float(np.mean(aris)) if aris else 1.0,
                })
            atac_pick = _rs.pick_stable_knee(atac_rows, tilt=atac_tilt, stability_floor=stability_floor)
            try:
                adata.close()
            except Exception:
                pass
        except Exception as e:
            log_event(run_dir, {"stage": "s7_clustering", "event": "atac_sweep_failed", "error": str(e)})
            atac_pick = {"resolution": 0.8, "rationale": "default fallback after sweep failure", "stable_region": []}

    # Persist sweep parquet
    sweep_df = pd.concat([
        pd.DataFrame(rna_rows).assign(modality="rna"),
        pd.DataFrame(atac_rows).assign(modality="atac") if atac_rows else pd.DataFrame(),
    ], ignore_index=True)
    sweep_df.to_parquet(art / "sweep.parquet")

    # Record recommended resolutions in provenance (source=derived; user can revise)
    if has_rna and rna_pick.get("resolution") is not None:
        _prov.set_param(params_path, "s7_clustering.rna.resolution",
                        float(rna_pick["resolution"]),
                        source="derived", confidence="medium",
                        rationale=rna_pick["rationale"],
                        method={"name": "resolution_sweep.stable_knee",
                                "code_ref": "executor/methods/resolution_sweep.py",
                                "inputs": {"grid": rna_grid, "seeds": seeds,
                                            "tilt": rna_tilt,
                                            "stability_floor": stability_floor}})
    if has_atac and atac_rows and atac_pick.get("resolution") is not None:
        _prov.set_param(params_path, "s7_clustering.atac.resolution",
                        float(atac_pick["resolution"]),
                        source="derived", confidence="medium",
                        rationale=atac_pick["rationale"],
                        method={"name": "resolution_sweep.stable_knee",
                                "code_ref": "executor/methods/resolution_sweep.py",
                                "inputs": {"grid": atac_grid, "seeds": seeds,
                                            "tilt": atac_tilt,
                                            "stability_floor": stability_floor}})

    # --- Adjacency comparison (policy-driven tie-breaker / interpretability check) ---
    adjacency_reports: dict[str, dict[str, Any]] = {}
    if has_rna and rna_rows:
        try:
            adjacency_reports["rna"] = _rcmp.adjacency_comparison(
                run_dir, modality="rna", grid=rna_grid,
                recommended=float(rna_pick["resolution"]), seed=seeds[0],
                sweep_rows=rna_rows,
            )
        except Exception as e:
            log_event(run_dir, {"stage": "s7_clustering", "event": "rna_adjacency_failed",
                                "error": str(e)})
    if has_atac and atac_rows:
        try:
            adjacency_reports["atac"] = _rcmp.adjacency_comparison(
                run_dir, modality="atac", grid=atac_grid,
                recommended=float(atac_pick["resolution"]), seed=seeds[0],
                sweep_rows=atac_rows,
            )
        except Exception as e:
            log_event(run_dir, {"stage": "s7_clustering", "event": "atac_adjacency_failed",
                                "error": str(e)})
    # Persist machine-readable adjacency report
    import json as _json
    (art / "adjacency_report.json").write_text(
        _json.dumps(adjacency_reports, indent=2, default=str)
    )

    # Human-readable recommendation + summary
    if branch == "paired":
        checkpoint_intro = (
            "# Clustering resolution review checkpoint\n"
            "\n"
            "**Paired multiome:** the sweep below is **diagnostic** — Leiden at each "
            "resolution is computed per modality on the shared joint cell set. Approved "
            "resolutions set `leiden_rna` / `leiden_atac` for UMAP colouring in the "
            "final `processed.h5mu`; they do **not** perform joint embedding or "
            "integrated clustering (WNN/MOFA+ are out of scope).\n"
        )
        approval_block = (
            "## Approval\n"
            "\n"
            "Confirm per-modality resolutions that yield interpretable diagnostic "
            "partitions (or revise), then approve before S8 runs:\n"
            "\n"
            "```bash\n"
            "Processing-MuAgent approve s7_clustering --config $CFG\n"
            "# or revise, e.g.:\n"
            "Processing-MuAgent revise s7_clustering s7_clustering.rna.resolution=1.2 --config $CFG\n"
            "```\n"
            "\n"
            "Final cluster labels are NOT assigned until approval.\n"
        )
    elif branch == "separate":
        checkpoint_intro = (
            "# Clustering resolution review checkpoint\n"
            "\n"
            "**Separate branch:** choose the RNA and ATAC resolutions that define "
            "**final** cluster labels in `rna_processed.h5ad` and `atac_processed.h5ad`. "
            "Use the sweep tables (n_clusters, silhouette, stability ARI) and optional "
            "adjacency comparison to decide.\n"
        )
        approval_block = (
            "## Approval\n"
            "\n"
            "Confirm the resolutions for your final processed outputs (or revise), "
            "then approve before S8 runs:\n"
            "\n"
            "```bash\n"
            "Processing-MuAgent approve s7_clustering --config $CFG\n"
            "```\n"
            "\n"
            "Final cluster labels are NOT assigned until approval.\n"
        )
    else:
        checkpoint_intro = (
            "# Clustering resolution review checkpoint\n"
            "\n"
            f"**{branch} branch:** choose the resolution that defines **final** cluster "
            "labels in the processed output. Use the sweep table below.\n"
        )
        approval_block = (
            "## Approval\n"
            "\n"
            "Confirm the resolution for your final processed output (or revise), "
            "then approve before S8 runs:\n"
            "\n"
            "```bash\n"
            "Processing-MuAgent approve s7_clustering --config $CFG\n"
            "```\n"
            "\n"
            "Final cluster labels are NOT assigned until approval.\n"
        )

    summary_lines = [
        checkpoint_intro,
        "",
        "## RNA",
    ]
    if has_rna and rna_rows:
        summary_lines += [
            f"- Recommended resolution: **{rna_pick['resolution']}**",
            f"- Rationale: {rna_pick['rationale']}",
            "",
            "| resolution | n_clusters | silhouette | stability ARI |",
            "|---:|---:|---:|---:|",
        ]
        for r in rna_rows:
            summary_lines.append(
                f"| {r['resolution']:.2f} | {r['n_clusters']} | "
                f"{r.get('silhouette', float('nan')):.3f} | {r['seed_stability_ari']:.3f} |"
            )
    else:
        summary_lines.append(f"- (not applicable — `{branch}` has no RNA modality)")
    summary_lines += ["", "## ATAC"]
    if has_atac and atac_rows:
        summary_lines.append(f"- Recommended resolution: **{atac_pick['resolution']}**")
        summary_lines.append(f"- Rationale: {atac_pick['rationale']}")
        summary_lines.append("")
        summary_lines.append("| resolution | n_clusters | silhouette | stability ARI |")
        summary_lines.append("|---:|---:|---:|---:|")
        for r in atac_rows:
            summary_lines.append(
                f"| {r['resolution']:.2f} | {r['n_clusters']} | "
                f"{r.get('silhouette', float('nan')):.3f} | {r['seed_stability_ari']:.3f} |"
            )
    elif has_atac:
        summary_lines.append("- (ATAC sweep failed — see log.jsonl)")
    else:
        summary_lines.append(f"- (not applicable — `{branch}` has no ATAC modality)")

    # Adjacency block (tie-breaker / interpretability check)
    summary_lines += ["", "## Adjacency comparison (tie-breaker / interpretability check)"]
    if "rna" in adjacency_reports:
        summary_lines.append(_rcmp.render_adjacency_block(adjacency_reports["rna"]))
    if "atac" in adjacency_reports:
        summary_lines.append(_rcmp.render_adjacency_block(adjacency_reports["atac"]))

    summary_lines += [approval_block]
    # Resolution summary is a user-facing deliverable — write directly to
    # deliverables/checkpoint/resolution_review/. The sweep.parquet + adjacency_report.json remain
    # under internal/artifacts/ as machine-readable intermediates.
    from ..run_paths import RunPaths
    summary_out = RunPaths(run_dir).resolution_summary_md
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.write_text("\n".join(summary_lines) + "\n")

    return {"rna_resolution": rna_pick["resolution"], "atac_resolution": atac_pick["resolution"]}


# ---------------------------------------------------------------------------
# Execute: read approved resolutions, assign final labels
# ---------------------------------------------------------------------------

def execute(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s7_clustering"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"
    branch = _prov.current_branch(str(params_path))
    has_rna = branch in ("paired", "separate", "rna_only")
    has_atac = branch in ("paired", "separate", "atac_only")

    p = plan["stages"]["s7_clustering"]["parameters"]
    seeds = list(p["seeds"]["value"])

    # Read approved (potentially user-revised) resolutions from parameters.yaml
    rna_res = _prov.get_value(params_path, "s7_clustering.rna.resolution")
    atac_res = _prov.get_value(params_path, "s7_clustering.atac.resolution")

    # --- RNA final labels ---
    if has_rna:
        if rna_res is None:
            raise RuntimeError("s7_clustering.rna.resolution not set in parameters.yaml; "
                                "did the propose step run?")
        rna = ad.read_h5ad(run_dir / "internal" / "artifacts" / "s6_dimred" / "rna_dimred.h5ad")
        sc.tl.leiden(rna, resolution=float(rna_res), random_state=seeds[0],
                     key_added="leiden_rna")
        _io.write_h5ad_safe(rna, art / "rna_clustered.h5ad")
    else:
        import scipy.sparse as sp
        _io.write_h5ad_safe(ad.AnnData(X=sp.csr_matrix((0, 0))), art / "rna_clustered.h5ad")

    # --- ATAC final labels ---
    if has_atac:
        if atac_res is None:
            raise RuntimeError("s7_clustering.atac.resolution not set in parameters.yaml; "
                                "did the propose step run?")
        try:
            import snapatac2 as snap
            atac_h5 = run_dir / "internal" / "artifacts" / "s5_atac_spectral" / "atac_spectral.h5ad"
            if not atac_h5.exists():
                atac_h5 = run_dir / "internal" / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"
            adata = snap.read(str(atac_h5))
            snap.tl.leiden(adata, resolution=float(atac_res), random_state=seeds[0],
                            key_added="leiden_atac")
            try:
                leiden_col = adata.obs["leiden_atac"]
            except Exception:
                leiden_col = np.asarray(adata.obs["leiden_atac"])
            pd.DataFrame({
                "barcode": [str(x) for x in adata.obs_names],
                "leiden_atac": np.asarray(leiden_col).astype(str),
            }).to_parquet(art / "atac_leiden_labels.parquet", index=False)
            try:
                adata.close()
            except Exception:
                pass
        except Exception as e:
            log_event(run_dir, {"stage": "s7_clustering", "event": "atac_finalize_failed",
                                "error": str(e)})
            raise

    log_event(run_dir, {"stage": "s7_clustering", "event": "done",
                         "rna_resolution": float(rna_res) if rna_res is not None else None,
                         "atac_resolution": float(atac_res) if atac_res is not None else None,
                         "source": "approved"})
    return {"rna_resolution": float(rna_res) if rna_res is not None else None,
            "atac_resolution": float(atac_res) if atac_res is not None else None}


# Backwards-compatible single-call entry (used by old Snakefile rule before split)
def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    propose(run_dir, plan)
    return execute(run_dir, plan)
