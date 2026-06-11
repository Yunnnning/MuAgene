"""Side-by-side Leiden resolution comparison helper.

Given an existing run directory with S6 PCA (RNA) + neighbor graph and S5 spectral ATAC, re-cluster
each modality at two user-specified resolutions and render a side-by-side UMAP
figure + summary table. Writes to <run_dir>/internal/artifacts/s7_clustering/
and surfaces comparison figures in deliverables/figures/.

Does not mutate the approved cluster labels on disk — this is a *comparison* tool
to inform resolution selection at the S7 checkpoint.

Also provides `adjacency_comparison()` — the policy-driven auto-comparison of the
recommended resolution against its nearest lower and higher neighbours on the sweep
grid. Outputs classifier verdicts (minor_reassignment / meaningful_split /
substantial_reorg / merge) and surface-condition flags (stability-tie,
silhouette-uninformative, small-medium-new-cluster).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from . import figures as _fig
from .atac_latent import ATAC_LATENT_KEY


FIGSIZE = (13, 5.5)

# ---------------------------------------------------------------------------
# Adjacency-comparison policy
# ---------------------------------------------------------------------------

# Size-threshold for "small/medium new cluster" detection (as fraction of total cells).
NEW_CLUSTER_MIN_FRAC = 0.05
NEW_CLUSTER_MAX_FRAC = 0.25

# ARI cut-offs for verdict classification.
ARI_MINOR = 0.90
ARI_MEANINGFUL = 0.75

# |Δ stability| cut-off for "adjacent stabilities are close" flag.
STABILITY_TIE_EPS = 0.05

# Silhouette range (max-min over grid) below which silhouette is considered uninformative.
SILHOUETTE_FLAT_RANGE = 0.02


def _cluster_sizes(labels: np.ndarray) -> dict[str, int]:
    s = pd.Series(labels.astype(str))
    return s.value_counts().sort_index().to_dict()


def _plot_two_panel(left_coords, left_labels, left_title,
                     right_coords, right_labels, right_title,
                     out_dir: Path, stem: str, suptitle: str) -> list[Path]:
    _fig._apply_style()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE)
    for ax, coords, labels, title in [
        (axes[0], left_coords, left_labels, left_title),
        (axes[1], right_coords, right_labels, right_title),
    ]:
        coords = np.asarray(coords)
        labels_arr = np.asarray(labels).astype(str)
        uniq = sorted(set(labels_arr), key=lambda s: (len(s), s))
        cmap = plt.get_cmap("tab20" if len(uniq) > 10 else "tab10")
        for i, v in enumerate(uniq):
            mask = labels_arr == v
            ax.scatter(coords[mask, 0], coords[mask, 1], s=8,
                       color=cmap(i % cmap.N), label=str(v),
                       alpha=0.85, linewidths=0)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(f"{title}  (n_clusters={len(uniq)})")
        ax.set_aspect("equal", adjustable="datalim")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(title="cluster", fontsize=_fig.FONT_SIZE - 3, markerscale=1.5,
                  bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0,
                  ncol=1 if len(uniq) <= 12 else 2)
    fig.suptitle(suptitle)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return _fig.save_figure(fig, out_dir, stem)


def compare_rna(run_dir: Path | str, resolutions: tuple[float, float], seed: int = 0) -> dict[str, Any]:
    run_dir = Path(run_dir)
    rna_path = run_dir / "internal" / "artifacts" / "s6_neighbors" / "rna_neighbors.h5ad"
    if not rna_path.exists():
        # fall back to the clustered output (also has X_umap + PCA + neighbors)
        rna_path = run_dir / "internal" / "artifacts" / "s7_clustering" / "rna_clustered.h5ad"
    rna = ad.read_h5ad(rna_path)
    # Ensure a UMAP is present
    if "X_umap" not in rna.obsm:
        sc.tl.umap(rna, random_state=seed)
    labels: dict[float, np.ndarray] = {}
    sizes: dict[float, dict[str, int]] = {}
    for res in resolutions:
        sc.tl.leiden(rna, resolution=float(res), random_state=seed,
                      key_added=f"_cmp_leiden_{res}")
        lbl = rna.obs[f"_cmp_leiden_{res}"].astype(str).to_numpy()
        labels[res] = lbl
        sizes[res] = _cluster_sizes(lbl)
    from .run_paths import RunPaths
    fig_out = RunPaths(run_dir).deliv_figures
    fig_out.mkdir(parents=True, exist_ok=True)
    figs = _plot_two_panel(
        rna.obsm["X_umap"], labels[resolutions[0]], f"RNA res={resolutions[0]}",
        rna.obsm["X_umap"], labels[resolutions[1]], f"RNA res={resolutions[1]}",
        out_dir=fig_out,
        stem=f"s7_compare_rna_res_{resolutions[0]}_vs_{resolutions[1]}",
        suptitle="RNA — Leiden resolution comparison",
    )
    return {"figures": [str(p) for p in figs],
            "sizes": {str(k): v for k, v in sizes.items()},
            "n_clusters": {str(k): len(v) for k, v in sizes.items()}}


def compare_atac(run_dir: Path | str, resolutions: tuple[float, float], seed: int = 0) -> dict[str, Any]:
    import snapatac2 as snap
    run_dir = Path(run_dir)
    atac_path = run_dir / "internal" / "artifacts" / "s5_atac_spectral" / "atac_spectral.h5ad"
    if not atac_path.exists():
        atac_path = run_dir / "internal" / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"
    adata = snap.read(str(atac_path))
    # Ensure ATAC UMAP exists
    try:
        coords = np.asarray(adata.obsm["X_umap"]) if "X_umap" in adata.obsm else None
    except Exception:
        coords = None
    if coords is None:
        try:
            snap.tl.umap(adata, random_state=seed, use_rep=ATAC_LATENT_KEY)
            coords = np.asarray(adata.obsm["X_umap"])
        except Exception:
            coords = None

    labels: dict[float, np.ndarray] = {}
    sizes: dict[float, dict[str, int]] = {}
    for res in resolutions:
        snap.tl.leiden(adata, resolution=float(res), random_state=seed,
                        key_added=f"_cmp_leiden_{res}")
        try:
            lbl = adata.obs[f"_cmp_leiden_{res}"].to_numpy().astype(str)
        except Exception:
            lbl = np.asarray(adata.obs[f"_cmp_leiden_{res}"]).astype(str)
        labels[res] = lbl
        sizes[res] = _cluster_sizes(lbl)
    try:
        adata.close()
    except Exception:
        pass
    if coords is None:
        return {"figures": [], "sizes": {str(k): v for k, v in sizes.items()},
                "n_clusters": {str(k): len(v) for k, v in sizes.items()},
                "note": "ATAC UMAP unavailable; no side-by-side figure rendered."}
    from .run_paths import RunPaths
    fig_out = RunPaths(run_dir).deliv_figures
    fig_out.mkdir(parents=True, exist_ok=True)
    figs = _plot_two_panel(
        coords, labels[resolutions[0]], f"ATAC res={resolutions[0]}",
        coords, labels[resolutions[1]], f"ATAC res={resolutions[1]}",
        out_dir=fig_out,
        stem=f"s7_compare_atac_res_{resolutions[0]}_vs_{resolutions[1]}",
        suptitle="ATAC — Leiden resolution comparison",
    )
    return {"figures": [str(p) for p in figs],
            "sizes": {str(k): v for k, v in sizes.items()},
            "n_clusters": {str(k): len(v) for k, v in sizes.items()}}


def render_comparison_summary(rna_report: dict, atac_report: dict) -> str:
    def sizes_line(sizes: dict[str, int]) -> str:
        items = sorted(sizes.items(), key=lambda kv: int(kv[0]) if kv[0].isdigit() else 10**9)
        return ", ".join(f"{c}:{n}" for c, n in items)

    out: list[str] = [
        "# Leiden resolution comparison",
        "",
        "_A side-by-side check to inform resolution selection at the S7 checkpoint._",
        "",
        "## RNA",
    ]
    for res, sizes in rna_report.get("sizes", {}).items():
        out.append(f"- **res={res}** → {rna_report['n_clusters'][res]} clusters | sizes: {sizes_line(sizes)}")
    out.append("")
    for p in rna_report.get("figures", []):
        out.append(f"- figure: `{Path(p).name}`")
    out.append("")
    out.append("## ATAC")
    for res, sizes in atac_report.get("sizes", {}).items():
        out.append(f"- **res={res}** → {atac_report['n_clusters'][res]} clusters | sizes: {sizes_line(sizes)}")
    if "note" in atac_report:
        out.append("")
        out.append(f"_{atac_report['note']}_")
    out.append("")
    for p in atac_report.get("figures", []):
        out.append(f"- figure: `{Path(p).name}`")
    out.append("")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# Adjacency-comparison (policy-driven auto-compare at S7 propose)
# ---------------------------------------------------------------------------


def _neighbours(grid: list[float], rec: float) -> tuple[float | None, float | None]:
    """Return (lower_neighbour, higher_neighbour) of rec in grid. Either may be None."""
    grid_s = sorted(set(grid))
    if rec not in grid_s:
        # Pick the closest grid value as the effective rec
        rec = min(grid_s, key=lambda g: abs(g - rec))
    idx = grid_s.index(rec)
    lower = grid_s[idx - 1] if idx > 0 else None
    higher = grid_s[idx + 1] if idx < len(grid_s) - 1 else None
    return lower, higher


def _verdict(n_a: int, n_b: int, ari: float, *, direction: str) -> str:
    """Classify the relationship between partitions A (lower/rec) and B (higher/rec).

    direction: "higher" if B is higher-resolution than A; "lower" otherwise.
    """
    if np.isnan(ari):
        return "substantial_reorg"
    if n_a == n_b:
        return "minor_reassignment" if ari >= ARI_MINOR else "substantial_reorg"
    if direction == "higher" and n_b > n_a:
        return "meaningful_split" if ari >= ARI_MEANINGFUL else "substantial_reorg"
    if direction == "lower" and n_b < n_a:
        return "merge" if ari >= ARI_MEANINGFUL else "substantial_reorg"
    return "substantial_reorg"


def _contingency(parent: np.ndarray, child: np.ndarray) -> dict[str, dict[str, int]]:
    """Cross-tab parent x child cluster membership counts."""
    p = pd.Series(np.asarray(parent).astype(str))
    c = pd.Series(np.asarray(child).astype(str))
    ct = pd.crosstab(p, c)
    return {str(row): {str(col): int(v) for col, v in ct.loc[row].items()} for row in ct.index}


def _new_clusters(parent: np.ndarray, child: np.ndarray) -> list[dict[str, Any]]:
    """Identify child clusters whose majority (>=75%) maps to a single parent cluster and
    whose size is in [NEW_CLUSTER_MIN_FRAC, NEW_CLUSTER_MAX_FRAC] of total cells.

    Returns a list of {label, size, frac, parent, parent_share}.
    """
    n_total = len(child)
    p = pd.Series(np.asarray(parent).astype(str))
    c = pd.Series(np.asarray(child).astype(str))
    out: list[dict[str, Any]] = []
    for child_label in sorted(c.unique(), key=lambda x: (len(x), x)):
        mask = c == child_label
        size = int(mask.sum())
        frac = size / n_total
        if frac < NEW_CLUSTER_MIN_FRAC or frac > NEW_CLUSTER_MAX_FRAC:
            continue
        parent_counts = p[mask].value_counts()
        if parent_counts.empty:
            continue
        top_parent = parent_counts.index[0]
        parent_share = parent_counts.iloc[0] / size
        if parent_share >= 0.75:
            out.append({
                "label": str(child_label),
                "size": size,
                "frac": round(frac, 4),
                "parent": str(top_parent),
                "parent_share": round(float(parent_share), 4),
            })
    return out


def _cluster_sizes_sorted(labels: np.ndarray) -> list[int]:
    s = pd.Series(np.asarray(labels).astype(str)).value_counts().sort_values(ascending=False)
    return [int(x) for x in s.values]


def _ari(a: np.ndarray, b: np.ndarray) -> float:
    from sklearn.metrics import adjusted_rand_score
    return float(adjusted_rand_score(a.astype(str), b.astype(str)))


def _rna_labels_at(adata, resolutions: list[float], seed: int) -> dict[float, np.ndarray]:
    out: dict[float, np.ndarray] = {}
    for res in resolutions:
        sc.tl.leiden(adata, resolution=float(res), random_state=seed,
                      key_added=f"_adj_leiden_{res}")
        out[res] = adata.obs[f"_adj_leiden_{res}"].astype(str).to_numpy()
    return out


def _atac_labels_at(adata, resolutions: list[float], seed: int) -> dict[float, np.ndarray]:
    import snapatac2 as snap
    out: dict[float, np.ndarray] = {}
    for res in resolutions:
        snap.tl.leiden(adata, resolution=float(res), random_state=seed,
                        key_added=f"_adj_leiden_{res}")
        try:
            lbl = adata.obs[f"_adj_leiden_{res}"].to_numpy()
        except Exception:
            lbl = np.asarray(adata.obs[f"_adj_leiden_{res}"])
        out[res] = np.asarray(lbl).astype(str)
    return out


def _compare_one(parent_labels: np.ndarray, child_labels: np.ndarray,
                  parent_res: float, child_res: float, direction: str) -> dict[str, Any]:
    n_parent = len(set(parent_labels))
    n_child = len(set(child_labels))
    ari = _ari(parent_labels, child_labels)
    verdict = _verdict(n_parent, n_child, ari, direction=direction)
    # "New" clusters only make sense when the child actually adds clusters (meaningful_split).
    # In minor_reassignment / substantial_reorg / merge, the 1:1 or many-to-one mapping is
    # not a genuine new-subpopulation signal.
    if direction == "higher" and verdict == "meaningful_split":
        new_clusters = _new_clusters(parent_labels, child_labels)
    else:
        new_clusters = []
    return {
        "parent_resolution": float(parent_res),
        "child_resolution": float(child_res),
        "direction": direction,
        "n_parent": n_parent,
        "n_child": n_child,
        "parent_sizes": _cluster_sizes_sorted(parent_labels),
        "child_sizes": _cluster_sizes_sorted(child_labels),
        "ari": round(ari, 4),
        "verdict": verdict,
        "new_clusters": new_clusters,
    }


def _surface_flags(sweep_rows: list[dict[str, Any]], rec: float,
                    lower: float | None, higher: float | None,
                    comparisons: list[dict[str, Any]]) -> list[str]:
    """Check surface-when conditions from the policy."""
    flags: list[str] = []
    by_res = {r["resolution"]: r for r in sweep_rows}
    # (a) adjacent stabilities close to recommendation
    stab_rec = by_res.get(rec, {}).get("seed_stability_ari", float("nan"))
    for neighbour in (lower, higher):
        if neighbour is None or neighbour not in by_res:
            continue
        stab_n = by_res[neighbour].get("seed_stability_ari", float("nan"))
        if not (np.isnan(stab_rec) or np.isnan(stab_n)):
            if abs(stab_rec - stab_n) < STABILITY_TIE_EPS:
                flags.append(f"stability_tie({neighbour} vs {rec}: "
                             f"Δ={abs(stab_rec - stab_n):.3f})")
    # (b) silhouette uninformative or missing
    sils = [r.get("silhouette", float("nan")) for r in sweep_rows]
    finite_sils = [s for s in sils if not (s is None or np.isnan(s))]
    if not finite_sils:
        flags.append("silhouette_missing")
    else:
        rng = max(finite_sils) - min(finite_sils)
        if rng < SILHOUETTE_FLAT_RANGE:
            flags.append(f"silhouette_uninformative(range={rng:.3f})")
    # (c) higher-res introduces small/medium new clusters
    for cmp in comparisons:
        if cmp["direction"] == "higher" and cmp["new_clusters"]:
            sizes = [nc["size"] for nc in cmp["new_clusters"]]
            flags.append(f"new_small_medium_clusters@res={cmp['child_resolution']}: "
                         f"{len(sizes)} cluster(s), sizes={sizes}")
    return flags


def adjacency_comparison(run_dir: Path | str, *, modality: str, grid: list[float],
                          recommended: float, seed: int,
                          sweep_rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Compare the recommended resolution to its nearest lower and higher grid values.

    Returns a dict with per-comparison verdicts, surface-condition flags, and (for
    convenience) the side-by-side figure paths. Writes figures under
    `deliverables/figures/` using the same naming convention as `compare_rna` /
    `compare_atac`.
    """
    run_dir = Path(run_dir)
    lower, higher = _neighbours(grid, recommended)
    comparisons: list[dict[str, Any]] = []
    figures: list[str] = []

    if modality == "rna":
        rna_path = run_dir / "internal" / "artifacts" / "s6_neighbors" / "rna_neighbors.h5ad"
        rna = ad.read_h5ad(rna_path)
        if "X_umap" not in rna.obsm:
            sc.tl.umap(rna, random_state=seed)
        needed = [r for r in (lower, recommended, higher) if r is not None]
        label_map = _rna_labels_at(rna, needed, seed=seed)
        coords = rna.obsm["X_umap"]
        umap_for_plot = coords
    elif modality == "atac":
        import snapatac2 as snap
        atac_path = run_dir / "internal" / "artifacts" / "s5_atac_spectral" / "atac_spectral.h5ad"
        if not atac_path.exists():
            atac_path = run_dir / "internal" / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"
        adata = snap.read(str(atac_path))
        # Ensure X_umap exists — S7 propose runs before S8 UMAP, so compute lazily.
        umap_for_plot = None
        try:
            if "X_umap" in adata.obsm:
                umap_for_plot = np.asarray(adata.obsm["X_umap"])
        except Exception:
            umap_for_plot = None
        if umap_for_plot is None:
            try:
                snap.tl.umap(adata, random_state=seed, use_rep=ATAC_LATENT_KEY)
                umap_for_plot = np.asarray(adata.obsm["X_umap"])
            except Exception:
                umap_for_plot = None
        needed = [r for r in (lower, recommended, higher) if r is not None]
        label_map = _atac_labels_at(adata, needed, seed=seed)
        try:
            adata.close()
        except Exception:
            pass
    else:
        raise ValueError(f"modality must be 'rna' or 'atac', got {modality!r}")

    # Adjacency check is numeric-only: ARI, cluster-count deltas, verdicts, and
    # surface flags. Per the output-layout policy, the adjacency check does NOT
    # save figures — the information is fully captured in `adjacency_report.json`
    # and the resolution-summary markdown block.
    rec_labels = label_map[recommended]
    if lower is not None:
        cmp_lo = _compare_one(label_map[lower], rec_labels,
                              parent_res=lower, child_res=recommended, direction="higher")
        # also produce the mirror "lower compared to rec" for merge detection:
        cmp_lo_mirror = _compare_one(rec_labels, label_map[lower],
                                      parent_res=recommended, child_res=lower, direction="lower")
        cmp_lo["merge_check"] = cmp_lo_mirror
        comparisons.append(cmp_lo)
    if higher is not None:
        cmp_hi = _compare_one(rec_labels, label_map[higher],
                              parent_res=recommended, child_res=higher, direction="higher")
        comparisons.append(cmp_hi)
    # `figures` intentionally remains an empty list — preserved in the schema for
    # downstream consumers that may still read the key.
    del umap_for_plot

    flags = _surface_flags(sweep_rows, recommended, lower, higher, comparisons)

    return {
        "modality": modality,
        "recommended": recommended,
        "lower": lower,
        "higher": higher,
        "comparisons": comparisons,
        "surface_flags": flags,
        "figures": figures,
        "policy": {
            "ari_minor": ARI_MINOR,
            "ari_meaningful": ARI_MEANINGFUL,
            "stability_tie_eps": STABILITY_TIE_EPS,
            "silhouette_flat_range": SILHOUETTE_FLAT_RANGE,
            "new_cluster_min_frac": NEW_CLUSTER_MIN_FRAC,
            "new_cluster_max_frac": NEW_CLUSTER_MAX_FRAC,
        },
    }


def render_adjacency_block(rep: dict[str, Any]) -> str:
    """Render the adjacency report as a markdown block suitable for resolution_summary.md."""
    modality = rep["modality"].upper()
    lines: list[str] = [
        f"### Adjacency check — {modality}",
        "",
        f"Recommended: **{rep['recommended']}**  (nearest lower: "
        f"{rep['lower']}, nearest higher: {rep['higher']})",
        "",
    ]
    if rep["surface_flags"]:
        lines.append("**Surface flags (this comparison is load-bearing for the decision):**")
        for f in rep["surface_flags"]:
            lines.append(f"- ⚠ {f}")
        lines.append("")
    for cmp in rep["comparisons"]:
        arrow = "→" if cmp["direction"] == "higher" else "←"
        lines.append(f"#### res={cmp['parent_resolution']}  {arrow}  res={cmp['child_resolution']}")
        lines.append(f"- n_clusters: {cmp['n_parent']} {arrow} {cmp['n_child']}")
        lines.append(f"- sorted sizes (parent): {cmp['parent_sizes']}")
        lines.append(f"- sorted sizes (child):  {cmp['child_sizes']}")
        lines.append(f"- ARI: {cmp['ari']}")
        lines.append(f"- **verdict: `{cmp['verdict']}`**")
        if cmp["new_clusters"]:
            lines.append("- new clusters (small/medium) in higher-res partition:")
            for nc in cmp["new_clusters"]:
                lines.append(f"  - label `{nc['label']}` — size {nc['size']} "
                             f"({nc['frac']*100:.1f}%), majority from parent "
                             f"`{nc['parent']}` ({nc['parent_share']*100:.0f}% share)")
        lines.append("")
    lines.append("_This is a tie-breaker / interpretability check, not a primary "
                 "selection metric. The selection still follows the stability-knee rule._")
    lines.append("")
    return "\n".join(lines)


def run_comparison(run_dir: Path | str,
                    rna_resolutions: tuple[float, float] = (1.0, 1.2),
                    atac_resolutions: tuple[float, float] = (0.6, 0.8)) -> Path:
    run_dir = Path(run_dir)
    rna_report = compare_rna(run_dir, rna_resolutions)
    atac_report = compare_atac(run_dir, atac_resolutions)
    out = run_dir / "internal" / "artifacts" / "s7_clustering" / "resolution_comparison.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render_comparison_summary(rna_report, atac_report))
    (out.parent / "resolution_comparison.json").write_text(json.dumps(
        {"rna": rna_report, "atac": atac_report,
         "rna_resolutions": list(rna_resolutions),
         "atac_resolutions": list(atac_resolutions)}, indent=2, default=str))
    return out
