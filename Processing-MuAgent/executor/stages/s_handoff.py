"""s_handoff — per-sample post-QC handoff bundle for Integration-MuAgent.

Emitted after post_qc_review approval (the split point between Preprocessing and
Integration). Writes two deliverables under deliverables/results/:

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

Integration-MuAgent reads the manifest to merge >=2 samples and re-count ATAC
fragments against a consensus peak set, so the prepared fragments + peaks BED must
survive the QC gate — see cli._cleanup_qc_intermediates (retain_for_integration).

Independently buildable terminal target (`run --target s_handoff`); orthogonal to
S4–S8, which still produce run_manifest.json. When S4–S8 move to Integration-MuAgent
this becomes Preprocessing's terminus.
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


def _load_mod(path: Path):
    """Load a post-doublet h5ad as an AnnData if it carries cells, else None.

    S3 always writes both rna_post_doublet.h5ad and atac_post_doublet.h5ad, but on
    single-modality branches the absent modality is a degenerate empty AnnData —
    those must not become a MuData modality.
    """
    if not path.exists():
        return None
    import anndata as ad
    try:
        adata = ad.read_h5ad(str(path))
    except Exception:
        return None
    if adata.n_obs == 0 or adata.n_vars == 0:
        return None
    return adata


def run(run_dir: Path | str, plan: dict[str, Any], workflow_branch: str) -> dict[str, Any]:
    """Assemble the per-sample post-QC .h5mu + manifest. `plan` is accepted for
    stage-signature parity (s3/s8) and currently unused; `workflow_branch` is the
    committed modality branch recorded by S0."""
    paths = RunPaths(Path(run_dir))
    run_dir = paths.run_dir
    params_path = str(paths.parameters_yaml)

    s3 = "s3_doublets"
    rna = _load_mod(paths.artifact(s3, "rna_post_doublet.h5ad"))
    atac = _load_mod(paths.artifact(s3, "atac_post_doublet.h5ad"))

    mods: dict[str, Any] = {}
    if rna is not None:
        mods["rna"] = rna
    if atac is not None:
        mods["atac"] = atac
    if not mods:
        raise RuntimeError(
            "s_handoff: no non-empty post-doublet modality under "
            f"{paths.artifact(s3, 'rna_post_doublet.h5ad').parent} — cannot assemble "
            "the post-QC handoff."
        )

    import mudata as mu
    mdata = mu.MuData(mods)
    h5mu_path = paths.deliv_results / f"post_qc_{run_dir.name}.h5mu"
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
    out = paths.deliv_results / "post_qc_manifest.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, default=str))

    log_event(run_dir, {"stage": "s_handoff", "event": "handoff_written",
                        "h5mu": rel(h5mu_path), "n_cells": n_cells})
    return {
        "post_qc_h5mu": rel(h5mu_path),
        "manifest": rel(out),
        "modality_branch": workflow_branch,
        "n_cells": n_cells,
    }
