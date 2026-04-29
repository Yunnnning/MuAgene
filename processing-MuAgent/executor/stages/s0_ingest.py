"""S0 — ingest, format detection, validation, pairing, metadata."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml

from .. import io as _io
from .. import pairing as _pair
from .. import metadata as _meta
from .. import provenance as _prov
from .. import hashing as _h
from ..log import log_event


def run(run_dir: Path | str, config: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    artifacts = run_dir / "internal" / "artifacts" / "s0_ingest"
    artifacts.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"

    rna_path = Path(config["rna_path"]) if config.get("rna_path") else None
    atac_frag_path = Path(config["atac_fragments_path"]) if config.get("atac_fragments_path") else None
    if rna_path is None and atac_frag_path is None:
        raise ValueError(
            "S0: at least one of `rna_path` or `atac_fragments_path` must be set in run.yaml."
        )

    # --- RNA side ---------------------------------------------------------
    rna = None
    rna_fmt: str | None = None
    single_file_multiome = False
    if rna_path is not None:
        rna_fmt = _io.detect_rna_format(rna_path)
        _prov.set_param(params_path, "ingest.rna_format", rna_fmt,
                        source="derived", confidence="high",
                        rationale=f"Autodetected from {rna_path}",
                        method={"name": "io.detect_rna_format",
                                "code_ref": "executor/io.py::detect_rna_format"})
        rna = _io.load_rna(rna_path, fmt=rna_fmt)
        if rna_fmt == "10x_h5":
            single_file_multiome = _io.detect_peaks_in_10x_h5(rna_path)

    # --- ATAC side --------------------------------------------------------
    frag_info: dict[str, Any] | None = None
    atac_bc: set[str] = set()
    genome_assembly = config.get("genome_assembly")
    if atac_frag_path is not None:
        if not genome_assembly:
            raise ValueError(
                "S0: `genome_assembly` is required for ATAC inputs; refusing to default."
            )
        frag_info = _io.validate_fragments(atac_frag_path)
        ok, msg = _io.cross_check_genome(set(frag_info["chromosomes"]), genome_assembly)
        if not ok:
            raise ValueError(f"S0 genome fingerprint mismatch: {msg}")
        atac_bc = _io.fragment_barcodes(atac_frag_path, limit=None)

    rna_bc: set[str] = set(rna.obs_names) if rna is not None else set()

    # Pairing — accepts empty sets on one side for single-modality branches.
    pr = _pair.detect_pairing(rna_bc, atac_bc, single_file_multiome=single_file_multiome)

    # --- Workflow branch derivation --------------------------------------
    if pr.status == "ambiguous":
        raise ValueError(
            f"S0 pairing is ambiguous (overlap={pr.overlap:.3f}); resolve before running preprocessing."
        )
    workflow_branch = pr.status  # paired | separate | rna_only | atac_only

    # User-declared branch (from `executor declare-branch`) — confirm or raise.
    declared = _prov.get_value(str(params_path), "plan.workflow_branch_declared", None)
    if declared is not None:
        if declared != workflow_branch:
            raise ValueError(
                f"S0: declared workflow_branch={declared!r} conflicts with detected "
                f"{workflow_branch!r}. Correct either the declaration or the inputs."
            )
        # Declaration matches detection — commit with source=user (schema forbids
        # method on source=user, so the confirmation is recorded in rationale).
        _prov.set_param(params_path, "plan.workflow_branch", workflow_branch,
                        source="user", confidence="high",
                        rationale=(f"User declared {declared!r} via `executor declare-branch`; "
                                   f"S0 detection via {pr.method} confirmed overlap={pr.overlap:.4f}."))
    else:
        _prov.set_param(params_path, "plan.workflow_branch", workflow_branch,
                        source="derived", confidence=pr.confidence,
                        rationale=f"From pairing status={pr.status}",
                        method={"name": "derive_workflow_branch",
                                "code_ref": "executor/stages/s0_ingest.py"})

    _prov.set_param(params_path, "ingest.pairing_decision",
                    {"status": pr.status, "confidence": pr.confidence, "method": pr.method,
                     "overlap": pr.overlap},
                    source="derived", confidence=pr.confidence,
                    rationale=f"Detected via {pr.method}; overlap={pr.overlap:.4f}",
                    method={"name": pr.method, "code_ref": "executor/pairing.py::detect_pairing"})

    # --- Metadata handling -----------------------------------------------
    meta_source = "reconstructed"
    meta_conf = "low"
    user_meta_df = None
    if config.get("metadata_path"):
        mpath = Path(config["metadata_path"])
        if mpath.exists():
            user_meta_df = _meta.load_user_metadata(mpath)
            join_col, coverage = _meta.identify_join_key(user_meta_df, rna_bc, atac_bc)
            meta_source = "provided"
            meta_conf = "high" if coverage >= 0.99 else "medium"
            _prov.set_param(params_path, "ingest.metadata_join_key", join_col,
                            source="derived", confidence=meta_conf,
                            rationale=f"Picked column with coverage {coverage:.3f}",
                            method={"name": "metadata.identify_join_key",
                                    "code_ref": "executor/metadata.py::identify_join_key"})

    _meta.reconstruct_minimal(rna_bc, atac_bc, artifacts / "metadata_minimal.tsv")

    _prov.set_param(params_path, "ingest.metadata_source", meta_source,
                    source="derived" if meta_source == "reconstructed" else "user",
                    confidence=meta_conf,
                    rationale="Minimal reconstruction from barcode union." if meta_source == "reconstructed" else "User-supplied metadata.",
                    method={"name": "metadata.minimal_reconstruction" if meta_source == "reconstructed" else "metadata.user_supplied",
                            "code_ref": "executor/metadata.py"} if meta_source == "reconstructed" else None)
    _prov.set_param(params_path, "ingest.metadata_unrecoverable",
                    _meta.unrecoverable_categories(meta_source),
                    source="derived", confidence="high",
                    rationale="Categories requiring user-supplied metadata.",
                    method={"name": "metadata.warn_unrecoverable",
                            "code_ref": "executor/metadata.py"})

    # --- Input hashes -> state.yaml --------------------------------------
    state_path = run_dir / "internal" / "state.yaml"
    state: dict[str, Any] = {}
    if state_path.exists():
        with state_path.open() as f:
            state = yaml.safe_load(f) or {}
    state.setdefault("input_hashes", {})
    if rna_path is not None:
        state["input_hashes"][str(rna_path)] = _h.sha256_file(rna_path)
    if atac_frag_path is not None:
        state["input_hashes"][str(atac_frag_path)] = _h.sha256_file(atac_frag_path)
    with state_path.open("w") as f:
        yaml.safe_dump(state, f)

    # --- Validation report -----------------------------------------------
    report: dict[str, Any] = {
        "workflow_branch": workflow_branch,
        "pairing": pr.as_dict(),
        "genome_assembly": genome_assembly,
        "metadata_source": meta_source,
        "single_file_multiome": single_file_multiome,
    }
    if rna is not None:
        report["rna_format"] = rna_fmt
        report["rna_n_cells"] = int(rna.n_obs)
        report["rna_n_genes"] = int(rna.n_vars)
    if frag_info is not None:
        report["atac_fragment_peek"] = frag_info
        report["atac_n_unique_barcodes"] = len(atac_bc)
    (artifacts / "validation_report.json").write_text(json.dumps(report, indent=2, default=str))

    # --- RNA ingest h5ad (always declared as an s0 output for DAG stability;
    #     written as an empty placeholder for atac_only so downstream rules'
    #     branch-aware input functions don't need to special-case existence). --
    rna_out = artifacts / "rna_ingest.h5ad"
    if rna is not None:
        if workflow_branch == "paired":
            common_bc = rna_bc & atac_bc
            rna_ingest = rna[rna.obs_names.isin(common_bc)].copy()
        else:
            rna_ingest = rna.copy()
        rna_ingest.layers["counts"] = rna_ingest.X.copy()
        rna_ingest.write_h5ad(rna_out)
    else:
        import scipy.sparse as sp
        import anndata as _ad
        _ad.AnnData(X=sp.csr_matrix((0, 0))).write_h5ad(rna_out)

    # --- ATAC ingest metadata (only if ATAC present) ---------------------
    if atac_frag_path is not None:
        (artifacts / "atac_ingest.json").write_text(json.dumps({
            "fragments_path": str(atac_frag_path),
            "tbi_path": str(Path(str(atac_frag_path) + ".tbi")),
            "barcodes_n": len(atac_bc),
            "chromosomes": sorted(set((frag_info or {}).get("chromosomes", []))),
        }, indent=2))

    log_event(run_dir, {"stage": "s0_ingest", "event": "done",
                        "workflow_branch": workflow_branch,
                        "n_cells_rna": int(rna.n_obs) if rna is not None else 0,
                        "n_barcodes_atac": len(atac_bc)})
    return report
