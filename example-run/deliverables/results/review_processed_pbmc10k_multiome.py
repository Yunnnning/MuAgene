"""Review processed_pbmc10k_multiome.h5mu — paired script version of the notebook.

Run from any working directory; RUN_DIR is baked in at generation time and
can be overridden via the PMA_RUN_DIR environment variable.
"""

# %% [markdown]
# # Review — `processed_pbmc10k_multiome.h5mu`
#
# Run directory: `.`

# %% [markdown]
# ## 1. Load data

# %%
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mudata as mu
import anndata as ad
from IPython.display import display

RUN_DIR = Path(os.environ.get("PMA_RUN_DIR", "."))
PROCESSED_DIR = RUN_DIR / "deliverables" / "results"
H5MU = PROCESSED_DIR / "processed_pbmc10k_multiome.h5mu"
RNA_H5AD = PROCESSED_DIR / "rna_processed.h5ad"
ATAC_H5AD = PROCESSED_DIR / "atac_processed.h5ad"

if H5MU.exists():
    mdata = mu.read_h5mu(H5MU)
elif RNA_H5AD.exists() and ATAC_H5AD.exists():
    # Unpaired branch: two independent h5ads
    mdata = mu.MuData({"rna": ad.read_h5ad(RNA_H5AD),
                       "atac": ad.read_h5ad(ATAC_H5AD)})
elif RNA_H5AD.exists():
    # rna_only branch
    mdata = mu.MuData({"rna": ad.read_h5ad(RNA_H5AD)})
elif ATAC_H5AD.exists():
    # atac_only branch
    mdata = mu.MuData({"atac": ad.read_h5ad(ATAC_H5AD)})
else:
    raise FileNotFoundError(
        f"No processed outputs found under {PROCESSED_DIR}")

print(mdata)

# %% [markdown]
# ## 2. Inspect contents

# %%
print(f"Modalities: {list(mdata.mod.keys())}")
for name, ad_ in mdata.mod.items():
    print(f"\n=== {name} ===")
    print(f"  shape:     {ad_.shape}")
    print(f"  obs cols:  {list(ad_.obs.columns)}")
    print(f"  obsm keys: {list(ad_.obsm.keys())}")
    if ad_.layers:
        print(f"  layers:    {list(ad_.layers.keys())}")

# Identify cluster-label columns
rna_cluster_col = next((c for c in mdata["rna"].obs.columns
                        if c.startswith("leiden_rna")), None)
atac_cluster_col = None
if "atac" in mdata.mod:
    atac_cluster_col = next((c for c in mdata["atac"].obs.columns
                             if c.startswith("leiden_atac")), None)

print(f"\nRNA cluster column:  {rna_cluster_col}")
print(f"ATAC cluster column: {atac_cluster_col}")

# %% [markdown]
# ### Cluster sizes

# %%
def _cluster_size_table(ad_, label_col: str) -> pd.DataFrame:
    counts = ad_.obs[label_col].value_counts().sort_index()
    return (counts.rename_axis(label_col)
                  .reset_index(name="n_cells"))

if rna_cluster_col:
    rna_sizes = _cluster_size_table(mdata["rna"], rna_cluster_col)
    print(f"RNA cluster sizes ({len(rna_sizes)} clusters, {int(rna_sizes['n_cells'].sum())} cells):")
    display(rna_sizes)

if atac_cluster_col:
    atac_sizes = _cluster_size_table(mdata["atac"], atac_cluster_col)
    print(f"ATAC cluster sizes ({len(atac_sizes)} clusters, {int(atac_sizes['n_cells'].sum())} cells):")
    display(atac_sizes)

# %% [markdown]
# ## 3. Clustering

# %% [markdown]
# ### UMAPs — coloured by cluster (reproduced from `processed_pbmc10k_multiome.h5mu`)

# %%
def _plot_umap_from_obj(ad_, coord_key: str, label_col: str, title: str) -> None:
    if coord_key not in ad_.obsm or label_col not in ad_.obs:
        print(f"(skip {title}: {coord_key} or {label_col} missing)")
        return
    coords = np.asarray(ad_.obsm[coord_key])
    labels = ad_.obs[label_col].astype(str)
    uniq = sorted(labels.unique(), key=lambda s: (len(s), s))
    cmap = plt.get_cmap("tab20" if len(uniq) > 10 else "tab10")
    fig, ax = plt.subplots(figsize=(6, 5))
    for i, lab in enumerate(uniq):
        mask = (labels == lab).values
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   s=8, color=cmap(i % cmap.N), label=lab,
                   alpha=0.85, linewidths=0)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left",
              fontsize=8, markerscale=2)
    fig.tight_layout()
    plt.show()

if rna_cluster_col:
    _plot_umap_from_obj(
        mdata["rna"],
        "X_umap_rna" if "X_umap_rna" in mdata["rna"].obsm else "X_umap",
        rna_cluster_col,
        "RNA UMAP — reproduced from processed_pbmc10k_multiome.h5mu",
    )
if atac_cluster_col:
    _plot_umap_from_obj(
        mdata["atac"],
        "X_umap_atac" if "X_umap_atac" in mdata["atac"].obsm else "X_umap",
        atac_cluster_col,
        "ATAC UMAP — reproduced from processed_pbmc10k_multiome.h5mu",
    )

# %% [markdown]
# ### Try a different resolution

# %%
# Try a DIFFERENT clustering resolution and regenerate the processed outputs.
#
# Re-running Leiden at a new resolution only changes the cluster *labels*; the
# UMAP embedding (computed from the latent representation — RNA: X_pca,
# ATAC: X_spectral) is unchanged. `recluster()` recomputes the labels, OVERWRITES
# the canonical `leiden_<modality>` column in memory, and recolours the UMAP.
# Nothing is written to disk until you call `save_processed()` — which regenerates
# the processed output file(s) so downstream tools use the new resolution.
import scanpy as sc

def recluster(modality="rna", resolution=0.7):
    if modality not in mdata.mod:
        print(f"(modality {modality!r} not present)"); return
    adata = mdata[modality]
    rep = "X_pca" if modality == "rna" else (
        "X_spectral" if "X_spectral" in adata.obsm else "X_lsi")
    if rep not in adata.obsm:
        print(f"(latent representation {rep!r} missing for {modality}; cannot re-cluster)")
        return
    sc.pp.neighbors(adata, use_rep=rep)
    label_col = f"leiden_{modality}"
    sc.tl.leiden(adata, resolution=resolution, key_added=label_col)
    adata.uns[f"{label_col}_resolution"] = float(resolution)
    n = adata.obs[label_col].nunique()
    print(f"{modality.upper()} re-clustered at resolution={resolution}: {n} clusters")
    print(f"  -> canonical column '{label_col}' updated in memory; "
          "call save_processed() to write it to disk.")
    umap_key = f"X_umap_{modality}" if f"X_umap_{modality}" in adata.obsm else "X_umap"
    _plot_umap_from_obj(adata, umap_key, label_col,
                        f"{modality.upper()} UMAP — Leiden res={resolution}")
    display(adata.obs[label_col].value_counts().sort_index()
            .rename_axis(label_col).reset_index(name="n_cells"))

def save_processed(*, overwrite=False, suffix="reclustered"):
    """Persist the re-clustered object(s) to disk so downstream tools use the new
    labels. Default: writes new file(s) next to the originals (non-destructive,
    `<name>_reclustered.<ext>`). Pass overwrite=True to replace the canonical
    processed output in place."""
    def _dest(orig):
        return orig if overwrite else orig.with_name(f"{orig.stem}_{suffix}{orig.suffix}")
    written = []
    if H5MU.exists():
        mdata.update()  # sync per-modality obs into the MuData before writing
        d = _dest(H5MU); mdata.write(str(d)); written.append(d)
    else:
        if RNA_H5AD.exists() and "rna" in mdata.mod:
            d = _dest(RNA_H5AD); mdata["rna"].write(str(d)); written.append(d)
        if ATAC_H5AD.exists() and "atac" in mdata.mod:
            d = _dest(ATAC_H5AD); mdata["atac"].write(str(d)); written.append(d)
    for w in written:
        print(f"Wrote {w}")
    return written

# Edit the resolution and re-run this cell. Examples:
recluster("rna", resolution=0.7)
# recluster("atac", resolution=0.5)
# Persist the re-clustered labels (non-destructive — writes *_reclustered files):
# save_processed()
# Or replace the canonical processed output in place:
# save_processed(overwrite=True)

