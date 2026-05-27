"""run_manifest.json writer — the handoff artifact."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import HANDOFF_CONTRACT_VERSION
from . import hashing as _h
from . import provenance as _prov
from .run_paths import RunPaths


def build_manifest(run_dir: Path | str, config: dict[str, Any]) -> dict[str, Any]:
    paths = RunPaths(Path(run_dir))
    run_dir = paths.run_dir

    def rel(p: str | Path) -> str:
        p = Path(p)
        try:
            return str(p.relative_to(run_dir))
        except ValueError:
            return str(p)

    def sha(p: Path) -> str | None:
        return _h.sha256_file(p) if p.exists() else None

    params_path = paths.parameters_yaml
    plan_path = paths.artifact("p2_plan", "preprocessing_plan.json")
    context_path = paths.artifact("p1_context", "context_extraction.json")

    # Outputs vary by branch — canonical locations are under deliverables/
    outputs: dict[str, Any] = {}
    if paths.processed_h5mu.exists():
        outputs["processed_h5mu"] = rel(paths.processed_h5mu)
    if paths.rna_processed_h5ad.exists():
        outputs["rna_processed_h5ad"] = rel(paths.rna_processed_h5ad)
    if paths.atac_processed_h5ad.exists():
        outputs["atac_processed_h5ad"] = rel(paths.atac_processed_h5ad)
    # User-facing figures (QC + UMAP) live in deliverables/figures/.
    figures = (
        [rel(f) for f in paths.deliv_figures.glob("*.png")]
        if paths.deliv_figures.exists() else []
    )
    outputs["figures"] = figures

    # Parameter hash
    param_hash = sha(params_path)

    # Input hashes from state.yaml if present
    input_hashes: dict[str, Any] = {}
    state_path = paths.state_yaml
    if state_path.exists():
        import yaml
        with state_path.open() as f:
            st = yaml.safe_load(f) or {}
        input_hashes = st.get("input_hashes", {})

    # Plan hash
    plan_hash = sha(plan_path)

    # Pairing decision (declared vs committed branch + method + downgrade reason).
    pairing_decision = _prov.get_value(str(params_path), "ingest.pairing_decision", {}) or {}
    pairing_block = {
        "declared": pairing_decision.get("declared"),
        "committed": pairing_decision.get("committed") or config.get("workflow_branch"),
        "method": pairing_decision.get("method"),
        "overlap": pairing_decision.get("overlap"),
        "reason": pairing_decision.get("downgrade_reason"),
    }
    # Final per-modality barcode counts — opened directly from the output files
    # so this records the actual state of what shipped, not what some upstream
    # stage claimed. n_joint is set only on the paired branch (h5mu).
    final_barcode_counts: dict[str, int | None] = {"rna": None, "atac": None, "joint": None}
    try:
        if paths.processed_h5mu.exists():
            import mudata as _mu
            mdata = _mu.read_h5mu(str(paths.processed_h5mu), backed="r")
            try:
                if "rna" in mdata.mod:
                    final_barcode_counts["rna"] = int(mdata.mod["rna"].n_obs)
                if "atac" in mdata.mod:
                    final_barcode_counts["atac"] = int(mdata.mod["atac"].n_obs)
                final_barcode_counts["joint"] = int(mdata.n_obs)
            finally:
                try:
                    mdata.file.close()
                except Exception:
                    pass
        else:
            import anndata as _ad
            if paths.rna_processed_h5ad.exists():
                rad = _ad.read_h5ad(str(paths.rna_processed_h5ad), backed="r")
                final_barcode_counts["rna"] = int(rad.n_obs)
                try:
                    rad.file.close()
                except Exception:
                    pass
            if paths.atac_processed_h5ad.exists():
                aad = _ad.read_h5ad(str(paths.atac_processed_h5ad), backed="r")
                final_barcode_counts["atac"] = int(aad.n_obs)
                try:
                    aad.file.close()
                except Exception:
                    pass
    except Exception:
        # Manifest must never block on count-reads; leave Nones if reads fail.
        pass

    manifest: dict[str, Any] = {
        "run_id": config.get("run_id"),
        "workflow_branch": config.get("workflow_branch"),
        "inputs": {
            "rna": {
                "path": config.get("rna_path"),
                "format": config.get("rna_format"),
                "sha256": input_hashes.get(config.get("rna_path")),
            },
            "atac_fragments": {
                "path": config.get("atac_fragments_path"),
                "sha256": input_hashes.get(config.get("atac_fragments_path")),
            },
            "metadata": {
                "path": config.get("metadata_path"),
                "source": config.get("metadata_source"),
            },
            "barcode_translation": {
                "path": config.get("barcode_translation_path"),
            },
            "atac_peaks": {
                "path": config.get("atac_peaks_path"),
            },
            "cell_metadata": {
                "path": config.get("cell_metadata_path"),
            },
        },
        "biological_context": {"ref": rel(context_path) if context_path.exists() else None},
        "preprocessing_plan": {
            "ref": rel(plan_path) if plan_path.exists() else None,
            "plan_hash": plan_hash,
        },
        "outputs": outputs,
        "pairing": pairing_block,
        "final_barcode_counts": final_barcode_counts,
        "parameters": {"ref": rel(params_path), "sha256": param_hash},
        "env": {
            "tool_versions": _h.tool_versions(),
            "seed": config.get("seed", 0),
        },
        "warnings": config.get("warnings", []),
        "handoff_contract_version": HANDOFF_CONTRACT_VERSION,
    }
    return manifest


def write_manifest(run_dir: Path | str, config: dict[str, Any]) -> Path:
    paths = RunPaths(Path(run_dir))
    m = build_manifest(run_dir, config)
    out = paths.run_manifest_json
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        json.dump(m, f, indent=2, default=str)
    return out
