"""Single-modality branches in `detect_pairing` — rna_only, atac_only, errors."""
from __future__ import annotations

import pytest

from executor.pairing import detect_pairing


def test_rna_only_empty_atac() -> None:
    r = detect_pairing({"AAA", "BBB", "CCC"}, set())
    assert r.status == "rna_only"
    assert r.confidence == "high"
    assert r.method == "pairing.rna_only_input"
    assert r.n_rna == 3 and r.n_atac == 0 and r.n_shared == 0
    assert r.overlap == 0.0


def test_atac_only_empty_rna() -> None:
    r = detect_pairing(set(), {"AAA", "BBB"})
    assert r.status == "atac_only"
    assert r.confidence == "high"
    assert r.method == "pairing.atac_only_input"
    assert r.n_rna == 0 and r.n_atac == 2 and r.n_shared == 0


def test_both_empty_raises() -> None:
    with pytest.raises(ValueError, match="both RNA and ATAC"):
        detect_pairing(set(), set())


def test_paired_exact_still_works() -> None:
    r = detect_pairing({"A", "B", "C"}, {"A", "B", "C"})
    assert r.status == "paired"
    assert r.method == "pairing.exact_barcode_match"


def test_separate_low_overlap_still_works() -> None:
    r = detect_pairing({"A", "B", "C"}, {"X", "Y", "Z"})
    assert r.status == "separate"
    assert r.method == "pairing.no_match"
