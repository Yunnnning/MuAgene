"""qc_handoff — per-sample post-QC handoff bundle for Integration-MuAgent.

Emitted after post_qc_review approval (the split point between Preprocessing and
Integration). Writes two deliverables under deliverables/qc/:

  post_qc_<run>.h5mu      — MuData{rna?, atac?} of POST-QC, POST-doublet cells,
                            UN-normalized (raw counts in rna.layers['counts']).
                            Integration does HVG/normalize with a batch_key, so the
                            handoff must NOT be pre-normalized. On single-modality
                            branches only the present mod is written.
  post_qc_manifest.json   — the cross-package handoff contract (schema
                            muagene.post_qc_handoff/1, versioned by
                            HANDOFF_CONTRACT_VERSION): modality branch, genome,
                            per-mod cell counts, and pointers to the RETAINED peaks
                            BED + prepared ATAC fragments.

ATAC encoding: S3 writes atac_post_doublet.h5ad with SnapATAC2's native, Blosc-
compressed writer, which plain anndata cannot decode without the HDF5 Blosc filter
plugin (absent in the standard env). Reading the ATAC side therefore goes through
snap.read (its Rust reader has the codec), and the mod is re-encoded as a *lean,
portable* AnnData — the raw fragments + chrom sizes + per-cell QC, gzip-compressed —
so the .h5mu is readable anywhere without the plugin and a downstream snap.read +
add_tile_matrix can rebuild the tile matrix. See _load_atac_mod.

Integration-MuAgent reads the manifest to merge >=2 samples and re-count ATAC
fragments against a consensus peak set, so the prepared fragments + peaks BED must
survive the QC gate — see cli._cleanup_qc_intermediates (retain_for_integration).

Independently buildable terminal target (`run --target qc_handoff`); orthogonal to
S4–S8, which still produce run_manifest.json. When S4–S8 move to Integration-MuAgent
this becomes Preprocessing's terminus. NOT a localrule — on HPC it runs as a SLURM
job because the ATAC materialisation is too heavy for the login/head node.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .. import HANDOFF_CONTRACT_VERSION
from .. import hashing as _h
from .. import io as _io
from .. import provenance as _prov
from ..log import log_event
from ..run_paths import RunPaths

# Prepared, chr-normalised fragment caches (either naming variant), in whichever
# stage dir they landed in. Retained past the gate for the consensus re-count.
_FRAG_NAMES = ("atac_fragments_cbf_chrnorm.tsv.gz", "atac_fragments_cbf.tsv.gz")
# Per-sample peak sets called at S2 (MACS3 preferred, ARC fallback).
_PEAK_NAMES = ("peaks_macs3.bed", "peaks_arc.bed")


def _find_first(paths: RunPaths, stages: tuple[str, ...], names: tuple[str, ...]) -> Path | None:
    for stage in stages:
        for name in names:
            p = paths.artifact(stage, name)
            if p.exists():
                return p
    return None


def _load_rna_mod(path: Path):
    """Load a post-doublet RNA h5ad as an AnnData if it carries cells, else None.

    S3 always writes rna_post_doublet.h5ad, but on atac_only the RNA side is a
    degenerate empty placeholder — that must not become a MuData modality. A genuine
    read error is re-raised (never silently swallowed into a partial handoff).
    """
    if not path.exists():
        return None
    import anndata as ad
    adata = ad.read_h5ad(str(path))
    if adata.n_obs == 0 or adata.n_vars == 0:
        return None
    return adata


def _load_atac_mod(path: Path):
    """Load the post-doublet ATAC modality for the handoff as a lean, portable AnnData.

    S3 writes atac_post_doublet.h5ad with SnapATAC2 (Blosc-compressed); read it via
    snap.read (Rust codec) and re-encode a portable AnnData carrying only what a
    downstream re-tile needs: the fragment insertion matrices (obsm, e.g.
    'fragment_paired'/'fragment_single'), the chrom sizes (uns['reference_sequences'],
    polars -> pandas), and the per-cell QC obs — with an EMPTY X. The derived,
    genome-scale tile matrix is intentionally not stored (snap.pp.add_tile_matrix
    rebuilds it from the fragments). This keeps the .h5mu small and plugin-free.

    Returns an anndata.AnnData, or None for the degenerate empty placeholder written
    on single-modality branches.
    """
    if not path.exists():
        return None
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    import snapatac2 as snap

    sdata = snap.read(str(path))  # backed; Rust reader handles the Blosc codec
    try:
        n_obs = int(sdata.n_obs)
        if n_obs == 0:
            return None
        obs = sdata.obs[:]
        obs_pd = obs.to_pandas() if hasattr(obs, "to_pandas") else pd.DataFrame(obs)
        obs_pd.index = [str(b) for b in list(sdata.obs_names)]
        adata = ad.AnnData(X=sp.csr_matrix((n_obs, 0), dtype=np.float32), obs=obs_pd)
        # Preserve every fragment/insertion obsm verbatim (paired- or single-end).
        for key in list(sdata.obsm.keys()):
            adata.obsm[key] = sdata.obsm[key]
        # reference_sequences (chrom sizes) is required by add_tile_matrix.
        rs = sdata.uns["reference_sequences"]
        adata.uns["reference_sequences"] = rs.to_pandas() if hasattr(rs, "to_pandas") else rs
        return adata
    finally:
        try:
            sdata.close()
        except Exception:
            pass


def run(run_dir: Path | str, plan: dict[str, Any], workflow_branch: str) -> dict[str, Any]:
    """Assemble the per-sample post-QC .h5mu + manifest. `plan` is accepted for
    stage-signature parity (s3/s8) and currently unused; `workflow_branch` is the
    committed modality branch recorded by S0."""
    paths = RunPaths(Path(run_dir))
    run_dir = paths.run_dir
    params_path = str(paths.parameters_yaml)

    s3 = "s3_doublets"
    has_rna = workflow_branch in ("paired", "separate", "rna_only")
    has_atac = workflow_branch in ("paired", "separate", "atac_only")

    rna = _load_rna_mod(paths.artifact(s3, "rna_post_doublet.h5ad")) if has_rna else None
    atac = _load_atac_mod(paths.artifact(s3, "atac_post_doublet.h5ad")) if has_atac else None

    # Loud failure: a modality the branch declares must be present. Never emit a
    # silent partial bundle (the old code swallowed an ATAC read error and dropped
    # the modality, producing an RNA-only h5mu on paired runs).
    if has_rna and rna is None:
        raise RuntimeError(
            f"qc_handoff: branch '{workflow_branch}' expects an RNA modality but "
            f"{paths.artifact(s3, 'rna_post_doublet.h5ad')} held no cells. Re-run S3 "
            "(e.g. revise an S1 threshold, which forces s1_rna_qc + s3_doublets to "
            "regenerate) before the handoff."
        )
    if has_atac and atac is None:
        raise RuntimeError(
            f"qc_handoff: branch '{workflow_branch}' expects an ATAC modality but "
            f"{paths.artifact(s3, 'atac_post_doublet.h5ad')} held no cells. Re-run S3 "
            "(e.g. revise an S2 threshold, which forces s2_atac_qc + s3_doublets to "
            "regenerate) before the handoff."
        )

    mods: dict[str, Any] = {}
    if rna is not None:
        mods["rna"] = rna
    if atac is not None:
        mods["atac"] = atac
    if not mods:
        raise RuntimeError(
            "qc_handoff: no non-empty post-doublet modality under "
            f"{paths.artifact(s3, 'rna_post_doublet.h5ad').parent} — cannot assemble "
            "the post-QC handoff."
        )

    import mudata as mu
    mdata = mu.MuData(mods)
    h5mu_path = paths.deliv_qc / f"post_qc_{run_dir.name}.h5mu"
    h5mu_path.parent.mkdir(parents=True, exist_ok=True)
    _io.write_mudata_safe(mdata, h5mu_path)

    # Retained ATAC artifacts the consensus re-count needs (kept past the gate via
    # retain_for_integration). May be None on RNA-only runs.
    peaks_bed = _find_first(paths, ("s2_atac_qc",), _PEAK_NAMES)
    fragments = _find_first(paths, ("qc_explore", "s2_atac_qc"), _FRAG_NAMES)
    add_chr_prefix = _prov.get_value(params_path, "s2_atac_qc.add_chr_prefix", None)
    genome_assembly = _prov.get_value(params_path, "ingest.genome_assembly", None)

    def rel(p: Path | None) -> str | None:
        if p is None:
            return None
        try:
            return str(Path(p).relative_to(run_dir))
        except ValueError:
            return str(p)

    n_cells = {
        "rna": int(rna.n_obs) if rna is not None else None,
        "atac": int(atac.n_obs) if atac is not None else None,
        "joint": int(mdata.n_obs),
    }

    manifest = {
        "schema": "muagene.post_qc_handoff/1",
        "handoff_contract_version": HANDOFF_CONTRACT_VERSION,
        "sample_run_dir": str(run_dir),
        "modality_branch": workflow_branch,
        "genome_assembly": genome_assembly,
        "post_qc_h5mu": rel(h5mu_path),
        "atac": {
            "peaks_bed": rel(peaks_bed),
            "fragments_prepared": rel(fragments),
            "add_chr_prefix": add_chr_prefix,
            # prepare_fragments_for_snapatac always emits UCSC-named fragments
            # (SnapATAC2 requirement), regardless of the source convention.
            "frag_chrom_convention": "ucsc" if fragments is not None else None,
        },
        "n_cells": n_cells,
        "parameters_ref": rel(paths.parameters_yaml),
        "tool_versions": _h.tool_versions(),
    }
    out = paths.deliv_qc / "post_qc_manifest.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, default=str))

    log_event(run_dir, {"stage": "qc_handoff", "event": "handoff_written",
                        "h5mu": rel(h5mu_path), "n_cells": n_cells})
    return {
        "post_qc_h5mu": rel(h5mu_path),
        "manifest": rel(out),
        "modality_branch": workflow_branch,
        "n_cells": n_cells,
    }
