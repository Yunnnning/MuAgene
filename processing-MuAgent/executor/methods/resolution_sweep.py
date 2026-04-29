"""Leiden resolution sweep with silhouette + seed-stability; stable-region knee picker."""
from __future__ import annotations

from typing import Any

import numpy as np
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, silhouette_score


def sweep_leiden(
    adata,
    *,
    grid: list[float],
    seeds: list[int],
    latent_key: str = "X_pca",
    neighbors_key: str = "neighbors",
) -> list[dict[str, Any]]:
    """Return per-resolution rows with n_clusters, silhouette, seed_stability."""
    rows: list[dict[str, Any]] = []
    latent = adata.obsm[latent_key]
    # If too many cells, subsample for silhouette to stay tractable
    n = adata.n_obs
    rng = np.random.default_rng(0)
    idx = rng.choice(n, size=min(n, 5000), replace=False)

    for res in grid:
        assignments = []
        for seed in seeds:
            sc.tl.leiden(adata, resolution=res, random_state=seed,
                         neighbors_key=neighbors_key, key_added=f"_tmp_leiden_{seed}")
            assignments.append(adata.obs[f"_tmp_leiden_{seed}"].astype(str).values)

        base = assignments[0]
        # silhouette on subsample
        try:
            unique = set(base[idx])
            if len(unique) > 1:
                sil = float(silhouette_score(latent[idx], base[idx], sample_size=None))
            else:
                sil = float("nan")
        except Exception:
            sil = float("nan")

        # seed stability: pairwise ARI mean
        aris: list[float] = []
        for i in range(len(seeds)):
            for j in range(i + 1, len(seeds)):
                aris.append(adjusted_rand_score(assignments[i], assignments[j]))
        stability = float(np.mean(aris)) if aris else 1.0

        rows.append({
            "resolution": float(res),
            "n_clusters": int(len(set(base))),
            "silhouette": sil,
            "seed_stability_ari": stability,
        })

    for seed in seeds:
        col = f"_tmp_leiden_{seed}"
        if col in adata.obs.columns:
            del adata.obs[col]
    return rows


def pick_stable_knee(rows: list[dict[str, Any]], *, tilt: str = "higher",
                     stability_floor: float = 0.85) -> dict[str, Any]:
    """Pick a resolution in the stable region.

    Stable region = contiguous set of resolutions with seed_stability_ari >= floor.
    Tilt 'higher' picks the highest resolution in the stable region (finer clusters);
    tilt 'lower' picks the lowest (coarser / avoids over-fragmentation).
    """
    if not rows:
        return {"resolution": None, "rationale": "empty sweep"}
    stable = [r for r in rows if r["seed_stability_ari"] >= stability_floor]
    if not stable:
        # fall back to best silhouette
        best = max(rows, key=lambda r: (-1.0 if np.isnan(r.get("silhouette") or float("nan")) else r["silhouette"]))
        return {
            "resolution": best["resolution"],
            "rationale": f"No resolution met stability floor {stability_floor}; fell back to max silhouette.",
            "stable_region": [],
        }
    stable_res = [r["resolution"] for r in stable]
    if tilt == "higher":
        chosen = max(stable_res)
    elif tilt == "lower":
        chosen = min(stable_res)
    else:
        chosen = sorted(stable_res)[len(stable_res) // 2]
    rationale = (
        f"Stable region {sorted(stable_res)} (ARI>={stability_floor}); "
        f"chose {chosen} per tilt={tilt}."
    )
    return {"resolution": chosen, "rationale": rationale, "stable_region": sorted(stable_res)}
