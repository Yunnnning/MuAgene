"""ATAC peak-export contracts.

Covers:
  1. `load_atac_feature_kind` reads S5's `feature_kind.txt` sidecar correctly
     for peak_matrix, tile_matrix, missing-file, and empty-file cases.
  2. Structural guards on s5_atac_lsi.py: peak-extraction branch is present,
     fallback to tile_matrix is present, and feature_kind.txt is written.
  3. Structural guard on io.py: `load_atac_from_10x_h5` exists (it's the
     single point s5 calls into for peak extraction).
"""
from __future__ import annotations

from pathlib import Path

import pytest

from executor.stages.s8_umap import load_atac_feature_kind


STAGE_S5 = Path("executor/stages/s5_atac_lsi.py")
STAGE_S8 = Path("executor/stages/s8_umap.py")
S5_RULE = Path("workflow/rules/s5_atac_lsi.smk")
IO_PY = Path("executor/io.py")


def _src(rel: Path) -> str:
    return (Path(__file__).parent.parent / rel).read_text()


# ---------------------------------------------------------------------------
# load_atac_feature_kind — sidecar reader
# ---------------------------------------------------------------------------


def test_kind_peak_matrix(tmp_path: Path) -> None:
    (tmp_path / "feature_kind.txt").write_text("peak_matrix\n")
    assert load_atac_feature_kind(tmp_path / "feature_kind.txt") == "peak_matrix"


def test_kind_tile_matrix(tmp_path: Path) -> None:
    (tmp_path / "feature_kind.txt").write_text("tile_matrix")
    assert load_atac_feature_kind(tmp_path / "feature_kind.txt") == "tile_matrix"


def test_kind_missing_file_returns_default(tmp_path: Path) -> None:
    # Default is an explicit unknown marker — never silently relabel a real
    # load as peak_matrix or tile_matrix when the sidecar is missing.
    assert load_atac_feature_kind(tmp_path / "missing.txt") == "unknown_feature_matrix"
    # Caller-supplied override still honored (e.g. for tests).
    assert load_atac_feature_kind(tmp_path / "missing.txt", default="latent_only") == "latent_only"


def test_kind_empty_file_returns_default(tmp_path: Path) -> None:
    (tmp_path / "feature_kind.txt").write_text("")
    assert load_atac_feature_kind(tmp_path / "feature_kind.txt") == "unknown_feature_matrix"


def test_kind_default_is_explicit_unknown() -> None:
    """Missing sidecars should not silently claim peak or tile semantics."""
    import inspect
    sig = inspect.signature(load_atac_feature_kind)
    default = sig.parameters["default"].default
    assert default not in {"tile_matrix", "peak_matrix"}, (
        f"load_atac_feature_kind default is {default!r}; it must not silently "
        "claim a concrete feature kind when the sidecar is missing."
    )
    assert default == "unknown_feature_matrix", (
        f"load_atac_feature_kind default is {default!r}; expected "
        "'unknown_feature_matrix' for missing/empty sidecars."
    )


# ---------------------------------------------------------------------------
# Structural contracts on S5
# ---------------------------------------------------------------------------


def test_s5_has_universal_peak_paths() -> None:
    """S5 must have both the ARC-h5 and MACS3 peak paths present."""
    src = _src(STAGE_S5)
    # Priority 1: ARC h5 path (gated on single_file_multiome from S0).
    assert "single_file_multiome" in src
    assert "io.load_atac_from_10x_h5" in src or "_io.load_atac_from_10x_h5" in src
    # Priority 2: MACS3 peak calling from fragments.
    assert "snap.tl.macs3" in src
    assert "snap.tl.merge_peaks" in src
    assert "snap.pp.make_peak_matrix" in src
    # S5 writes feature_kind.txt co-located with the matrix.
    assert "feature_kind.txt" in src


def test_s5_rule_proposal_mentions_peak_export_and_fallback() -> None:
    """Per-stage review text should mention the new export contract."""
    src = _src(S5_RULE)
    assert "spectral embedding" in src
    assert "peak matrix" in src
    assert "tile-matrix fallback" in src or "tile matrix fallback" in src


def test_s5_falls_back_to_tile_matrix_on_peak_failure() -> None:
    """Peak failure → explicit tile fallback path, NOT an interruption.

    Policy: preprocessing must not be interrupted by peak-export failures.
    S5 tries peak paths first, then falls back to the tile matrix that fed
    the spectral step, and honestly labels the fallback so downstream
    consumers can tell tile-fallback from a real peak export.
    """
    src = _src(STAGE_S5)
    # Fallback block is present and labelled.
    assert "tile-matrix fallback" in src.lower() or "tile_matrix fallback" in src.lower()
    # Fallback sets feature_kind="tile_matrix" (honest labelling).
    assert 'feature_kind = "tile_matrix"' in src
    # peak_source marker specifically for the fallback path.
    assert '"tile_matrix_fallback"' in src
    # An event is logged when the fallback engages (auditability).
    assert "tile_matrix_fallback_engaged" in src


def test_s5_does_not_raise_on_peak_failure() -> None:
    """The old hard-raise on peak failure must be gone; the run must continue."""
    src = _src(STAGE_S5)
    # The previous raise message content must not return — this is the exact
    # wording that tied a raise to peak-path failure. The phrase "peak
    # generation failed" is allowed in log_event messages (legit audit trail);
    # what must be absent is the old raise's specific sentence about export
    # being mandatory with no fallback.
    assert "ATAC feature-matrix export is mandatory" not in src, (
        "S5 must not assert peak export is mandatory; fallback is now explicit."
    )
    assert "no fallback to tile_matrix or latent_only is permitted" not in src, (
        "S5 must not carry the old 'no fallback permitted' wording; tile "
        "fallback is now the documented policy."
    )
    # Spectral/X_lsi is still a hard precondition (clustering can't proceed
    # without it) — that raise is unrelated to feature-matrix failures.
    assert '"X_lsi" not in adata.obsm' in src


def test_s5_feature_kind_matches_exported_source() -> None:
    """Exported kind must honestly reflect what was written.

    Three outcomes, three distinct labels:
      peak_matrix    — ARC or MACS3 peak path succeeded.
      tile_matrix    — tile-matrix fallback engaged.
      '' (empty)     — everything failed; S8 marks the output latent_only.
    """
    src = _src(STAGE_S5)
    # Both labels are set as literals in the code.
    assert 'feature_kind = "peak_matrix"' in src
    assert 'feature_kind = "tile_matrix"' in src
    # Empty sidecar path is also present (all-failed branch).
    assert '(art / "feature_kind.txt").write_text("")' in src
    # feature_kind is included in the lsi_summary JSON and in parameters.yaml.
    assert '"feature_kind": feature_kind' in src
    assert 's5_atac_lsi.feature_kind' in src


def test_s5_preserves_lsi_before_peak_export() -> None:
    """Tile/LSI clustering path is untouched; X_lsi must be written."""
    src = _src(STAGE_S5)
    # Spectral latent written to obsm.
    assert 'adata.obsm["X_lsi"]' in src
    # Post-export sanity: export refuses when X_lsi is missing.
    assert '"X_lsi" not in adata.obsm' in src


def test_s5_records_peak_source_provenance() -> None:
    """S5 records peak_source for real peak paths and for tile fallback."""
    src = _src(STAGE_S5)
    # peak_source is always recorded
    assert "s5_atac_lsi.peak_source" in src
    # Source values must appear in the code as literals.
    assert '"arc_h5"' in src
    assert '"macs3_from_fragments"' in src
    assert '"tile_matrix_fallback"' in src
    # peak_source_h5 is recorded conditionally (ARC path only).
    assert "s5_atac_lsi.peak_source_h5" in src


def test_s5_macs3_path_verifies_barcode_alignment() -> None:
    """MACS3 branch must check barcode identity + order, not just shape.

    Under-specification here would let `make_peak_matrix` silently drop or
    re-order cells and we'd export misaligned rows. The fix: verify every
    S5 barcode is present, reorder to S5's canonical order, then re-check.
    """
    src = _src(STAGE_S5)
    # Slice out the MACS3 block (between the Priority-2 and Priority-3 headers)
    # so we don't accidentally pick up ARC-path barcode checks or the tile
    # fallback.
    macs_marker = "Priority 2: MACS3 peak calling from fragments"
    tile_marker = "Priority 3: tile-matrix fallback"
    idx = src.find(macs_marker)
    end = src.find(tile_marker, idx)
    assert idx >= 0 and end > idx, "MACS3 block markers not found in s5 source"
    macs_block = src[idx:end]
    # Barcode existence check: MACS3 barcodes must cover every S5 barcode.
    assert "MACS3 peak matrix missing" in macs_block, (
        "MACS3 path must raise on missing barcodes, not rely on shape-only checks."
    )
    # Explicit reorder step.
    assert "order = [bc_to_idx[bc] for bc in s5_barcodes]" in macs_block, (
        "MACS3 path must reorder rows to s5_barcodes when identities match but order differs."
    )
    # Post-reorder identity check.
    assert "did not align to S5 barcodes after reorder" in macs_block, (
        "MACS3 path must verify final row identity after reordering."
    )


# ---------------------------------------------------------------------------
# Structural contract on io.py
# ---------------------------------------------------------------------------


def test_io_exposes_peak_loader() -> None:
    src = _src(IO_PY)
    assert "def load_atac_from_10x_h5" in src, (
        "io.load_atac_from_10x_h5 is the single entry point s5 uses for peak extraction."
    )
    # It must actually filter by feature type.
    assert "Peaks" in src
