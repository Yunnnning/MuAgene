import pytest
import anndata as ad
import numpy as np

from executor import pairing
from executor.stages import s0_ingest


def _pairing_result(status: str) -> pairing.PairingResult:
    return pairing.PairingResult(
        status=status,
        confidence="high",
        method="pairing.no_match",
        overlap=0.0,
        n_rna=10,
        n_atac=10,
        n_shared=0,
    )


def test_declared_paired_requires_confirmation_before_switching_to_unpaired():
    with pytest.raises(ValueError, match="explicitly declare `unpaired`"):
        s0_ingest._resolve_declared_branch("paired", _pairing_result("unpaired"))


def test_confirmed_pairing_keeps_paired_branch():
    assert s0_ingest._resolve_declared_branch(
        "paired", _pairing_result("paired")
    ) == "paired"


def test_disjoint_modalities_are_detected_as_unpaired():
    result = pairing.detect_pairing({"rna-1", "rna-2"}, {"atac-1", "atac-2"})
    assert result.status == "unpaired"


@pytest.mark.parametrize("sparse", [False, True])
def test_ingest_rejects_normalized_rna_counts(tmp_path, monkeypatch, sparse):
    from scipy import sparse as sp

    counts = np.array([[1.5, 2.0], [3.0, 4.25]])
    rna = ad.AnnData(sp.csr_matrix(counts) if sparse else counts)
    monkeypatch.setattr(s0_ingest._io, "detect_rna_format", lambda _: "h5ad")
    monkeypatch.setattr(s0_ingest._io, "detect_filtered_status", lambda *a, **k: "filtered")
    monkeypatch.setattr(s0_ingest._io, "load_rna_ingest", lambda *a, **k: (rna, None, None))

    with pytest.raises(ValueError, match="raw integer counts"):
        s0_ingest.run(tmp_path, {"rna_path": "normalized.h5ad"})
