"""S5 → S8 ATAC feature-matrix roundtrip + fallback contract.

Tests the loader used by S8 (`load_exported_atac_features`) directly, without
requiring SnapATAC2 / anndata objects. The loader is the single choke point for
S8's decision about whether to trust the S5 export; these tests pin its
validation rules:

  - Happy path: files present, aligned, names match columns → (mat, names).
  - Missing / invalid / misaligned inputs → (None, None), no silent repair.

Additional structural regression guards on the S5 writer + S8 fallback policy.
"""
from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pytest
import scipy.sparse as sp

from executor.stages.s8_umap import load_exported_atac_features


STAGE_S5 = Path("executor/stages/s5_atac_lsi.py")
STAGE_S8 = Path("executor/stages/s8_umap.py")


def _write_fixture(tmp_path: Path, *, mat, names: list[str]) -> tuple[Path, Path]:
    feat = tmp_path / "feature_matrix.npz"
    name_p = tmp_path / "feature_names.tsv"
    sp.save_npz(str(feat), mat)
    name_p.write_text("\n".join(names))
    return feat, name_p


def test_happy_path_valid_export(tmp_path: Path) -> None:
    mat = sp.csr_matrix(np.array([[1, 0, 2], [0, 3, 0], [4, 5, 0]], dtype=float))
    names = ["chr1:0-5000", "chr1:5000-10000", "chr2:0-5000"]
    feat, name_p = _write_fixture(tmp_path, mat=mat, names=names)
    out_mat, out_names = load_exported_atac_features(feat, name_p, n_obs_expected=3)
    assert out_mat is not None and out_names is not None
    assert out_mat.shape == (3, 3)
    assert out_names == names


def test_missing_files_returns_none(tmp_path: Path) -> None:
    assert load_exported_atac_features(
        tmp_path / "missing.npz", tmp_path / "missing.tsv", n_obs_expected=5,
    ) == (None, None)


def test_row_mismatch_returns_none(tmp_path: Path) -> None:
    mat = sp.csr_matrix(np.eye(3))
    feat, name_p = _write_fixture(tmp_path, mat=mat, names=["a", "b", "c"])
    assert load_exported_atac_features(feat, name_p, n_obs_expected=4) == (None, None)


def test_column_name_count_mismatch_returns_none(tmp_path: Path) -> None:
    mat = sp.csr_matrix(np.eye(3))
    # 2 names for 3 columns — must fail (no silent trim).
    feat, name_p = _write_fixture(tmp_path, mat=mat, names=["a", "b"])
    assert load_exported_atac_features(feat, name_p, n_obs_expected=3) == (None, None)


def test_zero_columns_returns_none(tmp_path: Path) -> None:
    # Degenerate export: 3 cells × 0 features. Not a usable representation.
    mat = sp.csr_matrix((3, 0))
    sp.save_npz(str(tmp_path / "feature_matrix.npz"), mat)
    (tmp_path / "feature_names.tsv").write_text("")
    assert load_exported_atac_features(
        tmp_path / "feature_matrix.npz", tmp_path / "feature_names.tsv",
        n_obs_expected=3,
    ) == (None, None)


def test_corrupt_npz_returns_none(tmp_path: Path) -> None:
    (tmp_path / "feature_matrix.npz").write_text("not a real npz")
    (tmp_path / "feature_names.tsv").write_text("a\nb\nc")
    assert load_exported_atac_features(
        tmp_path / "feature_matrix.npz", tmp_path / "feature_names.tsv",
        n_obs_expected=3,
    ) == (None, None)


# ---------------------------------------------------------------------------
# Structural guards on the S5 writer + S8 fallback policy — these pin the
# contract without running SnapATAC2.
# ---------------------------------------------------------------------------


def _src(path: Path) -> str:
    full = Path(__file__).parent.parent / path
    return full.read_text()


def test_s5_verifies_before_export() -> None:
    """S5 must verify every peak path (shape + names + X_lsi) before exporting."""
    src = _src(STAGE_S5)
    # Refuses to export without a computed spectral latent (clustering guarantee).
    assert re.search(r'"X_lsi" not in adata\.obsm', src)
    # Both peak paths validate shape against n_obs_expected and refuse on invalid shape.
    assert "invalid for n_obs=" in src
    # Both peak paths validate var_names count against matrix columns.
    assert "var_names length" in src
    # Shape × n_obs check — "shape" text appears in both validators.
    assert re.search(r"peak (matrix|var_names)", src, re.IGNORECASE)


def test_s8_fallback_is_zero_column_not_one_column() -> None:
    """S8 must never fabricate a 1-column .X placeholder."""
    src = _src(STAGE_S8)
    # The documented fallback shape: (n_obs, 0).
    assert re.search(r"sp\.csr_matrix\(\(len\(atac_barcodes\), 0\)\)", src), (
        "S8 fallback must use zero columns; a 1-column placeholder is misleading."
    )
    # The old misleading pattern must not return.
    assert not re.search(r"sp\.csr_matrix\(\(len\(atac_barcodes\), 1\)\)", src), (
        "S8 still uses a 1-column fake .X — remove it in favour of zero-column."
    )


def test_s8_feature_kind_is_gated_on_real_load() -> None:
    """S8 reads the kind from S5's sidecar only when a real matrix loaded; latent_only otherwise."""
    src = _src(STAGE_S8)
    # Real-load branch defers to S5's feature_kind.txt sidecar (peak_matrix or tile_matrix).
    assert 'load_atac_feature_kind(kind_path)' in src, (
        "S8 must read the kind from S5's sidecar, not hard-code a single value."
    )
    # Latent-only branch marks latent_only plainly.
    assert 'atac_adata.uns["atac_feature_kind"] = "latent_only"' in src
