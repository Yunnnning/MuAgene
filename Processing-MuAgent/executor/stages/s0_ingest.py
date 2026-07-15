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
from .. import translation as _translation
from ..log import log_event


def _resolve_declared_branch(
    declared: str | None,
    pairing_result: _pair.PairingResult,
) -> str:
    """Require explicit consent before changing a declared paired workflow."""
    if declared == "paired" and pairing_result.status != "paired":
        raise ValueError(
            "S0 could not validate the declared `paired` workflow "
            f"(detected {pairing_result.status!r}, overlap={pairing_result.overlap:.4f}). "
            "Review the pairing diagnostics, then provide a barcode translation, correct "
            "the inputs, or explicitly declare `unpaired` before rerunning S0."
        )
    return pairing_result.status


def run(run_dir: Path | str, config: dict[str, Any]) -> dict[str, Any]:
    run_dir = Path(run_dir)
    artifacts = run_dir / "internal" / "artifacts" / "s0_ingest"
    artifacts.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"

    rna_path = Path(config["rna_path"]) if config.get("rna_path") else None
    rna_raw_path = Path(config["rna_raw_path"]) if config.get("rna_raw_path") else None
    atac_frag_path = Path(config["atac_fragments_path"]) if config.get("atac_fragments_path") else None
    if rna_path is None and rna_raw_path is None and atac_frag_path is None:
        raise ValueError(
            "S0: at least one of `rna_path`, `rna_raw_path`, or `atac_fragments_path` "
            "must be set in run.yaml."
        )

    # --- RNA side ---------------------------------------------------------
    rna = None                 # cell-called / filtered matrix (used downstream)
    rna_raw_full = None        # full barcode set (only set when a raw matrix is
                               # supplied or detected) — required by SoupX.
    rna_fmt: str | None = None
    rna_filtered_status: str | None = None  # "filtered" | "raw"
    cell_calling_diag: dict[str, Any] | None = None
    single_file_multiome = False

    # If only `rna_raw_path` is supplied (no filtered matrix), promote it to
    # `rna_path` for downstream code: we'll cell-call internally.
    rna_input_path = rna_path or rna_raw_path
    if rna_input_path is not None:
        rna_fmt = _io.detect_rna_format(rna_input_path)
        rna_filtered_status = _io.detect_filtered_status(rna_input_path, fmt=rna_fmt)
        # If user explicitly used `rna_raw_path` to point at a filtered file,
        # honour the user's intent: treat as raw so cell calling runs.
        # (Detection threshold may misclassify small / unusual matrices.)
        if rna_path is None and rna_raw_path is not None:
            rna_filtered_status = "raw"
        _prov.set_param(params_path, "ingest.rna_format", rna_fmt,
                        source="derived", confidence="high",
                        rationale=f"Autodetected from {rna_input_path}",
                        method={"name": "io.detect_rna_format",
                                "code_ref": "executor/io.py::detect_rna_format"})
        _prov.set_param(params_path, "ingest.rna_filtered_status", rna_filtered_status,
                        source="derived", confidence="high",
                        rationale=("Detected via barcode count "
                                    f"(threshold={_io.RAW_BARCODE_THRESHOLD})."),
                        method={"name": "io.detect_filtered_status",
                                "code_ref": "executor/io.py::detect_filtered_status"})

        if rna_fmt == "10x_h5":
            single_file_multiome = _io.detect_peaks_in_10x_h5(rna_input_path)

        # Ingest the RNA matrix via the shared SSOT: load -> (raw input) deterministic
        # barcode-rank knee cell-calling -> add a `counts` layer. S1a reconstructs
        # rna_ingest.h5ad with the identical call after the post-QC cleanup deletes it,
        # so the two stages never diverge on what "ingested RNA" means.
        rna, rna_raw_full, cell_calling_diag = _io.load_rna_ingest(
            rna_input_path, fmt=rna_fmt, filtered_status=rna_filtered_status)

        if cell_calling_diag is not None:
            # Raw matrix → record the barcode-rank knee cell-calling provenance.
            _prov.set_param(params_path, "ingest.cell_calling_method",
                            cell_calling_diag.get("method", "unknown"),
                            source="derived", confidence="high",
                            rationale=("Raw RNA matrix supplied; cells called via "
                                        "barcode-rank knee on log-counts vs log-rank curve."),
                            method={"name": "io.call_cells_from_raw",
                                    "code_ref": "executor/io.py::call_cells_from_raw"})
            _prov.set_param(params_path, "ingest.cell_calling_threshold",
                            float(cell_calling_diag.get("threshold", 0.0)),
                            source="derived", confidence="medium",
                            rationale=("Total-count threshold derived from the "
                                        "barcode-rank knee."),
                            method={"name": "io.barcode_rank_knee",
                                    "code_ref": "executor/io.py::barcode_rank_knee"})
            _prov.set_param(params_path, "ingest.n_cells_called",
                            int(cell_calling_diag.get("n_kept", rna.n_obs)),
                            source="derived", confidence="high",
                            rationale="Cells retained after barcode-rank knee call.",
                            method={"name": "io.call_cells_from_raw",
                                    "code_ref": "executor/io.py::call_cells_from_raw"})
        elif rna_raw_path is not None:
            # Filtered input + optional companion raw matrix (filtered + raw both
            # supplied): used downstream by SoupX for soup-profile estimation.
            raw_fmt = _io.detect_rna_format(rna_raw_path)
            rna_raw_full = _io.load_rna(rna_raw_path, fmt=raw_fmt)
            _prov.set_param(params_path, "ingest.rna_raw_format", raw_fmt,
                            source="derived", confidence="high",
                            rationale=f"Autodetected from {rna_raw_path}",
                            method={"name": "io.detect_rna_format",
                                    "code_ref": "executor/io.py::detect_rna_format"})

        # Integer-counts guard. seurat_v3 HVG (S4) and Scrublet (S3) require
        # raw integer counts. Refuse early if X is already normalized.
        try:
            x = rna.X
            sample = x[:50] if hasattr(x, "shape") and x.shape[0] > 50 else x
            arr = sample.toarray() if hasattr(sample, "toarray") else np.asarray(sample)
            arr = arr.ravel()
            arr = arr[np.isfinite(arr)]
            if arr.size and not np.allclose(arr, np.round(arr), atol=1e-6):
                raise ValueError(
                    "S0: RNA matrix .X does not look like raw integer counts (sample "
                    "contains non-integer values). seurat_v3 HVG and Scrublet require "
                    "raw counts. If you supplied a normalized .h5ad, replace it with "
                    "a raw-counts matrix or move the raw counts into .X."
                )
        except ValueError:
            raise
        except Exception:
            pass

    # --- ATAC side --------------------------------------------------------
    frag_info: dict[str, Any] | None = None
    atac_bc: set[str] = set()
    atac_fragments_file_barcodes_n: int | None = None
    atac_n_cells_report = 0
    atac_cell_whitelist_path: Path | None = None
    genome_assembly = config.get("genome_assembly")
    if atac_frag_path is not None:
        if not genome_assembly:
            raise ValueError(
                "S0: `genome_assembly` is required for ATAC inputs; refusing to default."
            )
        atac_fmt = _io.detect_atac_format(atac_frag_path)
        if atac_fmt == "bed4":
            # Auto-convert 4-column BED to standard 5-column bgzipped fragments.
            # The derived file is written next to the source; the source is unchanged.
            atac_frag_path = _io.convert_bed4_to_fragments(atac_frag_path)
            _prov.set_param(params_path, "ingest.atac_format_original", "bed4",
                            source="derived", confidence="high",
                            rationale="4-column BED detected; converted to 5-column fragments.tsv.gz",
                            method={"name": "io.convert_bed4_to_fragments",
                                    "code_ref": "executor/io.py::convert_bed4_to_fragments"})
            _prov.set_param(params_path, "ingest.atac_fragments_derived_path",
                            str(atac_frag_path),
                            source="derived", confidence="high",
                            rationale="Path to converted fragments.tsv.gz used downstream.",
                            method={"name": "io.convert_bed4_to_fragments",
                                    "code_ref": "executor/io.py::convert_bed4_to_fragments"})
        frag_info = _io.validate_fragments(atac_frag_path)
        ok, msg = _io.cross_check_genome(set(frag_info["chromosomes"]), genome_assembly)
        if not ok:
            raise ValueError(f"S0 genome fingerprint mismatch: {msg}")
        _prov.set_param(params_path, "ingest.genome_assembly", genome_assembly,
                        source="user", confidence="high",
                        rationale="Declared in run.yaml; cross-checked against ATAC fragment chromosomes in S0.")
        atac_bc = _io.fragment_barcodes(atac_frag_path, limit=None)
        atac_fragments_file_barcodes_n = len(atac_bc)

    rna_bc: set[str] = set(rna.obs_names) if rna is not None else set()

    # Cell-level ATAC count for reporting. Cell Ranger ARC fragments.tsv.gz carries
    # GEX barcodes for every droplet; the filtered RNA matrix is the cell-called
    # subset. Do not restrict atac_bc here — pairing and downstream import still
    # see the full fragments barcode set.
    atac_n_cells_report = len(atac_bc)
    if rna_bc and atac_bc:
        rel_report = _pair.detect_subset_relation(rna_bc, atac_bc)
        if rel_report.relation == "rna_subset_of_atac":
            atac_n_cells_report = len(rna_bc)

    # When the fragments file carries all droplets but RNA defines the cell-called
    # subset (typical Cell Ranger ARC), persist an explicit whitelist for S2's
    # SnapATAC2 import so n_cells_pre matches the filtered matrix, not every
    # barcode with ATAC reads.
    if (
        rna_bc
        and atac_frag_path is not None
        and atac_fragments_file_barcodes_n is not None
        and atac_fragments_file_barcodes_n > atac_n_cells_report
    ):
        atac_cell_whitelist_path = artifacts / "atac_cell_barcodes.tsv"
        _io.write_text_safe(
            atac_cell_whitelist_path,
            "\n".join(sorted(rna_bc)) + "\n",
        )

    # Pairing — accepts empty sets on one side for single-modality branches.
    pr_initial = _pair.detect_pairing(rna_bc, atac_bc, single_file_multiome=single_file_multiome)

    # --- Diagnostics ladder for the workflow_branch decision -------------
    # Detection is advisory; declaration + supplied mappings drive the committed branch.
    # Ladder rungs (first hit wins; thresholds in executor/pairing.py):
    #   1. Direct Jaccard >= PAIRING_OVERLAP_THRESHOLD -> pairing.exact_barcode_match
    #   2. ATAC or RNA barcode subset (>= SUBSET_COVERAGE_THRESHOLD of smaller set)
    #   3. Suffix-normalized Jaccard or subset -> prefix_suffix_normalized / *_subset_of_*
    #   4. `barcode_translation_path` or cell_metadata translation table
    #   5. otherwise -> require explicit confirmation before switching to unpaired
    declared = _prov.get_value(str(params_path), "plan.workflow_branch_declared", None)
    barcode_translation_path = config.get("barcode_translation_path")
    cell_metadata_path = config.get("cell_metadata_path")

    committed_branch = pr_initial.status
    pairing_result = pr_initial
    def _ladder_step(pr: _pair.PairingResult, **extra: Any) -> dict[str, Any]:
        step: dict[str, Any] = {
            "step": pr.method,
            "status": pr.status,
            "overlap": pr.overlap,
        }
        if pr.subset_relation:
            step["subset_relation"] = pr.subset_relation
            step["subset_coverage"] = pr.subset_coverage
        step.update(extra)
        return step

    ladder_steps: list[dict[str, Any]] = [_ladder_step(pr_initial)]
    translation_table_loaded: dict[str, str] | None = None
    translation_source: str | None = None

    if declared == "paired" and committed_branch != "paired" and rna_bc and atac_bc:
        # Try barcode_translation_path first.
        if barcode_translation_path:
            tpath = Path(barcode_translation_path)
            if tpath.exists():
                try:
                    translation_table_loaded = _translation.load_translation_tsv(tpath)
                    translation_source = f"barcode_translation_path={tpath}"
                except Exception as e:
                    log_event(run_dir, {"stage": "s0_ingest",
                                         "event": "translation_table_load_failed",
                                         "path": str(tpath), "error": str(e)})

        # Then try cell_metadata_path if it carries both rna_barcode and atac_barcode.
        if translation_table_loaded is None and cell_metadata_path:
            mpath = Path(cell_metadata_path)
            if mpath.exists():
                try:
                    import pandas as _pd
                    head_df = _pd.read_csv(mpath, sep="\t", dtype=str, nrows=5,
                                           keep_default_na=False)
                    cols_l = {c.lower() for c in head_df.columns}
                    if "atac_barcode" in cols_l and (
                        "rna_barcode" in cols_l or "gex_barcode" in cols_l
                    ):
                        translation_table_loaded = _translation.load_translation_tsv(mpath)
                        translation_source = f"cell_metadata_path={mpath}"
                except Exception as e:
                    log_event(run_dir, {"stage": "s0_ingest",
                                         "event": "cell_metadata_translation_load_failed",
                                         "path": str(mpath), "error": str(e)})

        if translation_table_loaded:
            pr_trans = _pair.pairing_via_translation(
                rna_bc, atac_bc, translation_table_loaded,
            )
            ladder_steps.append(_ladder_step(
                pr_trans, source=translation_source,
                n_pairs=len(translation_table_loaded),
            ))
            if pr_trans.status == "paired":
                committed_branch = "paired"
                pairing_result = pr_trans
                _translation.write_translation_parquet(
                    translation_table_loaded,
                    artifacts / "barcode_translation.parquet",
                )
            else:
                committed_branch = "unpaired"
                pairing_result = pr_trans
        else:
            committed_branch = "unpaired"
            pairing_result = pr_initial

    committed_branch = _resolve_declared_branch(declared, pairing_result)

    # Ambiguous overlap with no resolution path is still a hard stop — needs human input.
    if pairing_result.status == "ambiguous" and committed_branch == "ambiguous":
        raise ValueError(
            f"S0 pairing is ambiguous (overlap={pairing_result.overlap:.3f}); resolve "
            "before running preprocessing — supply `barcode_translation_path` or declare "
            "the branch explicitly with `executor declare-branch`."
        )

    # Single-modality declarations that contradict the detected modality set are still
    # a hard error: rna_only with ATAC inputs (or atac_only with RNA inputs) signals
    # a data-hygiene problem the user must resolve.
    if declared is not None and declared != "paired" and declared != committed_branch:
        raise ValueError(
            f"S0: declared workflow_branch={declared!r} conflicts with detected "
            f"{committed_branch!r}. For the paired<->unpaired decision, supply a "
            "`barcode_translation_path`; for rna_only/atac_only declarations, either "
            "remove the unwanted modality from the inputs or correct the declaration."
        )

    workflow_branch = committed_branch  # paired | unpaired | rna_only | atac_only

    # Commit workflow_branch with appropriate provenance source/method.
    if declared == "paired" and workflow_branch == "paired" \
            and pairing_result.method == "pairing.translation_table":
        _prov.set_param(params_path, "plan.workflow_branch", workflow_branch,
                        source="user", confidence="high",
                        rationale=(f"User declared 'paired'; pairing established via "
                                   f"translation table ({translation_source}), "
                                   f"overlap={pairing_result.overlap:.4f}."))
    elif declared == "paired" and workflow_branch == "paired":
        _prov.set_param(params_path, "plan.workflow_branch", workflow_branch,
                        source="user", confidence="high",
                        rationale=(f"User declared 'paired' via `executor declare-branch`; "
                                   f"S0 detection via {pairing_result.method} confirmed "
                                   f"overlap={pairing_result.overlap:.4f}."))
    elif declared is not None and declared == workflow_branch:
        _prov.set_param(params_path, "plan.workflow_branch", workflow_branch,
                        source="user", confidence="high",
                        rationale=(f"User declared {declared!r} via `executor declare-branch`; "
                                   f"S0 detection via {pairing_result.method} matched."))
    else:
        _prov.set_param(params_path, "plan.workflow_branch", workflow_branch,
                        source="derived", confidence=pairing_result.confidence,
                        rationale=f"From pairing status={pairing_result.status}",
                        method={"name": "derive_workflow_branch",
                                "code_ref": "executor/stages/s0_ingest.py"})

    pairing_record: dict[str, Any] = {
        "status": pairing_result.status,
        "confidence": pairing_result.confidence,
        "method": pairing_result.method,
        "overlap": pairing_result.overlap,
        "ladder": ladder_steps,
        "declared": declared,
        "committed": workflow_branch,
        "thresholds": {
            "jaccard_paired": _pair.PAIRING_OVERLAP_THRESHOLD,
            "subset_coverage": _pair.SUBSET_COVERAGE_THRESHOLD,
            "ambiguous_low": _pair.AMBIGUOUS_OVERLAP_LOW,
        },
    }
    if pairing_result.subset_relation:
        pairing_record["subset_relation"] = pairing_result.subset_relation
        pairing_record["subset_coverage"] = pairing_result.subset_coverage
    if translation_source:
        pairing_record["translation_source"] = translation_source
        pairing_record["n_translation_pairs"] = len(translation_table_loaded or {})

    _prov.set_param(params_path, "ingest.pairing_decision", pairing_record,
                    source="derived", confidence=pairing_result.confidence,
                    rationale=f"Detected via {pairing_result.method}; overlap={pairing_result.overlap:.4f}",
                    method={"name": pairing_result.method,
                            "code_ref": "executor/pairing.py"})

    # Record the new optional input paths so manifest + reproducibility surface them.
    for cfg_key, param_key in (
        ("barcode_translation_path", "ingest.barcode_translation_path"),
        ("atac_peaks_path", "ingest.atac_peaks_path"),
        ("cell_metadata_path", "ingest.cell_metadata_path"),
    ):
        val = config.get(cfg_key)
        if val:
            _prov.set_param(params_path, param_key, str(val),
                            source="user", confidence="high",
                            rationale=f"Supplied by user in run.yaml as {cfg_key}.")

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
    if rna_raw_path is not None:
        state["input_hashes"][str(rna_raw_path)] = _h.sha256_file(rna_raw_path)
    if atac_frag_path is not None:
        state["input_hashes"][str(atac_frag_path)] = _h.sha256_file(atac_frag_path)
    with state_path.open("w") as f:
        yaml.safe_dump(state, f)

    # --- Validation report -----------------------------------------------
    report: dict[str, Any] = {
        "workflow_branch": workflow_branch,
        "pairing": pairing_record,
        "genome_assembly": genome_assembly,
        "metadata_source": meta_source,
        "single_file_multiome": single_file_multiome,
        "rna_filtered_status": rna_filtered_status,
        "has_raw_matrix": rna_raw_full is not None,
    }
    # Surface the optional input paths in the report so users see them
    # without having to grep parameters.yaml.
    for k in ("barcode_translation_path", "atac_peaks_path", "cell_metadata_path",
              "rna_path", "rna_raw_path"):
        if config.get(k):
            report[k] = str(config[k])
    if rna is not None:
        report["rna_format"] = rna_fmt
        report["rna_n_cells"] = int(rna.n_obs)
        report["rna_n_genes"] = int(rna.n_vars)
    if cell_calling_diag is not None:
        report["cell_calling"] = cell_calling_diag
    if rna_raw_full is not None:
        report["rna_raw_n_barcodes"] = int(rna_raw_full.n_obs)
    if frag_info is not None:
        report["atac_fragment_peek"] = frag_info
        report["atac_n_unique_barcodes"] = atac_n_cells_report
        if (
            atac_fragments_file_barcodes_n is not None
            and atac_fragments_file_barcodes_n != atac_n_cells_report
        ):
            report["atac_fragments_file_barcodes_n"] = atac_fragments_file_barcodes_n
            report["atac_barcodes_source"] = "rna_cell_call"
    _io.write_text_safe(artifacts / "validation_report.json", json.dumps(report, indent=2, default=str))

    # --- RNA ingest h5ad (always declared as an s0 output for DAG stability;
    #     written as an empty placeholder for atac_only so downstream rules'
    #     branch-aware input functions don't need to special-case existence). --
    rna_out = artifacts / "rna_ingest.h5ad"
    if rna is not None:
        # No pre-intersection at S0 — modality-specific QC (S1/S2) runs on each
        # modality's full barcode set. For the paired branch, S3 enforces the
        # joint barcode intersection after doublet removal. The `counts` layer was
        # already added by io.load_rna_ingest (the SSOT), so write `rna` directly —
        # this is a deletable cache (S1a reconstructs it via io.load_rna_ingest).
        _io.write_h5ad_safe(rna, rna_out)
    else:
        import scipy.sparse as sp
        import anndata as _ad
        _io.write_h5ad_safe(_ad.AnnData(X=sp.csr_matrix((0, 0))), rna_out)

    # --- Raw matrix input ref (symlink; used by SoupX in S1a) ------------
    # Never copy the full raw matrix into the run dir — point at the user path.
    raw_source_path: Path | None = None
    raw_source_fmt: str | None = None
    if rna_raw_path is not None:
        raw_source_path = rna_raw_path
        raw_source_fmt = _io.detect_rna_format(rna_raw_path)
    elif rna_filtered_status == "raw" and rna_input_path is not None:
        raw_source_path = rna_input_path
        raw_source_fmt = rna_fmt
    if raw_source_path is not None and raw_source_fmt is not None:
        _io.write_input_ref(artifacts / "rna_raw", raw_source_path, fmt=raw_source_fmt)
        _prov.set_param(params_path, "ingest.rna_raw_source_path", str(raw_source_path.resolve()),
                        source="derived", confidence="high",
                        rationale="Symlinked at S0; raw matrix is read from the original input path.",
                        method={"name": "io.write_input_ref",
                                "code_ref": "executor/io.py::write_input_ref"})
        legacy_raw = artifacts / "rna_raw.h5ad"
        if legacy_raw.exists():
            legacy_raw.unlink()

    # --- ATAC ingest metadata (only if ATAC present) ---------------------
    if atac_frag_path is not None:
        atac_ingest_meta: dict[str, Any] = {
            "fragments_path": str(atac_frag_path),
            "tbi_path": str(Path(str(atac_frag_path) + ".tbi")),
            "barcodes_n": atac_n_cells_report,
            "chromosomes": sorted(set((frag_info or {}).get("chromosomes", []))),
        }
        if (
            atac_fragments_file_barcodes_n is not None
            and atac_fragments_file_barcodes_n != atac_n_cells_report
        ):
            atac_ingest_meta["fragments_file_barcodes_n"] = atac_fragments_file_barcodes_n
            atac_ingest_meta["barcodes_source"] = "rna_cell_call"
        if atac_cell_whitelist_path is not None:
            atac_ingest_meta["cell_barcode_whitelist"] = str(atac_cell_whitelist_path)
        _io.write_text_safe(artifacts / "atac_ingest.json", json.dumps(atac_ingest_meta, indent=2))

    # --- Preprocessing plan assembly (merged in from the former p2_plan rule) ---
    # assemble_plan is deterministic and needs no heavy data — only the context,
    # the just-written ingest report, and the committed branch. Assembling it here
    # lets the single S0 job also run the QC exploration on the in-memory matrices.
    from .. import plan_assembler as _pa
    ctx_path = (run_dir / "internal" / "artifacts" / "p1_context"
                / "context_extraction.json")
    try:
        ctx = json.loads(ctx_path.read_text()) if ctx_path.exists() else {"fields": {}}
    except Exception:
        ctx = {"fields": {}}
    sample_type = (ctx.get("fields", {}).get("sample_type") or {}).get("value", "unknown")
    plan = _pa.assemble_plan(
        run_dir,
        workflow_branch=workflow_branch,
        sample_type=sample_type,
        ingest=report,
        s1a_ambient_method=config.get("s1a_ambient_method"),
    )
    _, plan_hash = _pa.write_plan(run_dir, plan)
    _prov.set_param(
        params_path, "plan.plan_hash", plan_hash,
        source="derived", confidence="high",
        rationale="sha256 of preprocessing_plan.json",
        method={"name": "sha256_bytes", "code_ref": "executor/hashing.py::sha256_bytes"},
    )

    # --- QC exploration on the in-memory data (no reload) ---
    # Pass the already-loaded RNA matrix so qc_explore never re-reads rna_ingest.h5ad;
    # ATAC fragments are imported once here. Best-effort — exploration failure must
    # not block plan review (a degraded report is recoverable).
    from .. import qc_explore as _qc_explore
    try:
        _qc_explore.run(run_dir, rna_adata=rna)
    except Exception as e:
        log_event(run_dir, {"stage": "s0_ingest", "event": "qc_explore_failed",
                            "error": str(e)})

    log_event(run_dir, {"stage": "s0_ingest", "event": "done",
                        "workflow_branch": workflow_branch,
                        "n_cells_rna": int(rna.n_obs) if rna is not None else 0,
                        "n_barcodes_atac": atac_n_cells_report,
                        "n_fragments_file_barcodes": atac_fragments_file_barcodes_n})
    return report
