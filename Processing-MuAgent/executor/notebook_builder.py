"""Generate a user-facing Jupyter review notebook for each run.

Writes:
    <run_dir>/deliverables/results/review_processed_<run>.ipynb
    <run_dir>/deliverables/results/review_processed_<run>.py

The notebook is self-contained: RUN_DIR is baked in at generation time and can be
overridden at runtime via the `PMA_RUN_DIR` environment variable so the notebook
remains portable if the run directory is relocated.

Design goals:
    - Concise, user-friendly, no debugging noise.
    - Three sections only: load data, inspect contents, clustering.
    - Reproduce the per-cluster UMAP from the loaded object (QC / ambient figures
      are NOT reproduced here — they already live in deliverables/figures/ and the
      QC-review checkpoint).
    - Let the user re-cluster at a different resolution and regenerate the
      processed output (new labels + recoloured UMAP).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Cell helpers — build nbformat v4.5 JSON directly (no nbformat dependency)
# ---------------------------------------------------------------------------

def _cell_id(source: str, salt: str) -> str:
    """Stable 12-char id derived from source + salt. nbformat ≥4.5 requires cell ids."""
    h = hashlib.sha1((salt + "::" + source).encode()).hexdigest()
    return h[:12]


def _md(source: str, salt: str = "") -> dict[str, Any]:
    return {
        "cell_type": "markdown",
        "id": _cell_id(source, salt or f"md:{source[:20]}"),
        "metadata": {},
        "source": source.splitlines(keepends=True),
    }


def _code(source: str, salt: str = "") -> dict[str, Any]:
    return {
        "cell_type": "code",
        "id": _cell_id(source, salt or f"code:{source[:20]}"),
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source.splitlines(keepends=True),
    }


# ---------------------------------------------------------------------------
# Cell content — kept small and composable for the paired .py export below
# ---------------------------------------------------------------------------

def _cell_header(run_dir: str, run_name: str) -> dict:
    h5mu_file = f"processed_{run_name}.h5mu"
    return _md(
        f"# Review — `{h5mu_file}`\n"
        f"\n"
        f"Run directory: `{run_dir}`\n"
    )


def _cell_load(run_name: str) -> str:
    h5mu_file = f"processed_{run_name}.h5mu"
    return f"""\
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mudata as mu
import anndata as ad
from IPython.display import display

RUN_DIR = Path(os.environ.get("PMA_RUN_DIR", "__BAKED_RUN_DIR__"))
PROCESSED_DIR = RUN_DIR / "deliverables" / "results"
H5MU = PROCESSED_DIR / "{h5mu_file}"
RNA_H5AD = PROCESSED_DIR / "rna_processed.h5ad"
ATAC_H5AD = PROCESSED_DIR / "atac_processed.h5ad"

if H5MU.exists():
    mdata = mu.read_h5mu(H5MU)
elif RNA_H5AD.exists() and ATAC_H5AD.exists():
    # Separate branch: two independent h5ads
    mdata = mu.MuData({{"rna": ad.read_h5ad(RNA_H5AD),
                       "atac": ad.read_h5ad(ATAC_H5AD)}})
elif RNA_H5AD.exists():
    # rna_only branch
    mdata = mu.MuData({{"rna": ad.read_h5ad(RNA_H5AD)}})
elif ATAC_H5AD.exists():
    # atac_only branch
    mdata = mu.MuData({{"atac": ad.read_h5ad(ATAC_H5AD)}})
else:
    raise FileNotFoundError(
        f"No processed outputs found under {{PROCESSED_DIR}}")

print(mdata)
"""


_CELL_INSPECT = """\
print(f"Modalities: {list(mdata.mod.keys())}")
for name, ad_ in mdata.mod.items():
    print(f"\\n=== {name} ===")
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

print(f"\\nRNA cluster column:  {rna_cluster_col}")
print(f"ATAC cluster column: {atac_cluster_col}")
"""


_CELL_CLUSTER_SIZES = """\
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
"""


def _cell_umap_repro(run_name: str) -> str:
    h5mu_file = f"processed_{run_name}.h5mu"
    return f"""\
def _plot_umap_from_obj(ad_, coord_key: str, label_col: str, title: str) -> None:
    if coord_key not in ad_.obsm or label_col not in ad_.obs:
        print(f"(skip {{title}}: {{coord_key}} or {{label_col}} missing)")
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
        "RNA UMAP — reproduced from {h5mu_file}",
    )
if atac_cluster_col:
    _plot_umap_from_obj(
        mdata["atac"],
        "X_umap_atac" if "X_umap_atac" in mdata["atac"].obsm else "X_umap",
        atac_cluster_col,
        "ATAC UMAP — reproduced from {h5mu_file}",
    )
"""


def _cell_recluster() -> str:
    return """\
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
    \"\"\"Persist the re-clustered object(s) to disk so downstream tools use the new
    labels. Default: writes new file(s) next to the originals (non-destructive,
    `<name>_reclustered.<ext>`). Pass overwrite=True to replace the canonical
    processed output in place.\"\"\"
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
"""


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def _cells(run_dir: str, run_name: str) -> list[dict[str, Any]]:
    h5mu_file = f"processed_{run_name}.h5mu"

    def bake(s: str) -> str:
        return s.replace("__BAKED_RUN_DIR__", run_dir)

    return [
        _cell_header(run_dir, run_name),
        _md("## 1. Load data"),
        _code(bake(_cell_load(run_name))),
        _md("## 2. Inspect contents"),
        _code(_CELL_INSPECT),
        _md("### Cluster sizes"),
        _code(_CELL_CLUSTER_SIZES),
        _md("## 3. Clustering"),
        _md(f"### UMAPs — coloured by cluster (reproduced from `{h5mu_file}`)"),
        _code(_cell_umap_repro(run_name)),
        _md("### Try a different resolution"),
        _code(_cell_recluster()),
    ]


def build_notebook(run_dir: Path | str) -> dict[str, Any]:
    run_dir = Path(run_dir)
    return {
        "cells": _cells(str(run_dir), run_dir.name),
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "name": "python",
                "mimetype": "text/x-python",
                "file_extension": ".py",
                "pygments_lexer": "ipython3",
                "nbconvert_exporter": "python",
                "version": "3.10",
                "codemirror_mode": {"name": "ipython", "version": 3},
            },
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def build_script(run_dir: Path | str) -> str:
    """Paired .py version with jupytext percent-format cell markers."""
    run_dir = Path(run_dir)
    h5mu_file = f"processed_{run_dir.name}.h5mu"
    header = (
        f'"""Review {h5mu_file} — paired script version of the notebook.\n'
        "\n"
        "Run from any working directory; RUN_DIR is baked in at generation time and\n"
        "can be overridden via the PMA_RUN_DIR environment variable.\n"
        '"""\n\n'
    )
    parts: list[str] = [header]
    for cell in _cells(str(run_dir), run_dir.name):
        src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
        if cell["cell_type"] == "markdown":
            parts.append("# %% [markdown]\n")
            for line in src.splitlines():
                parts.append(f"# {line}\n" if line else "#\n")
            parts.append("\n")
        else:
            parts.append("# %%\n")
            parts.append(src if src.endswith("\n") else src + "\n")
            parts.append("\n")
    return "".join(parts)


def write_review_notebook(run_dir: Path | str) -> tuple[Path, Path]:
    from .run_paths import RunPaths
    run_dir = Path(run_dir).resolve()
    paths = RunPaths(run_dir)
    paths.deliv_results.mkdir(parents=True, exist_ok=True)

    ipynb = paths.review_notebook_ipynb
    ipynb.write_text(json.dumps(build_notebook(run_dir), indent=1))

    py = paths.review_notebook_py
    py.write_text(build_script(run_dir))

    return ipynb, py
