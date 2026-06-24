"""When rna_ingest.h5ad has been cleaned (post-QC cleanup), S1a reconstructs it
deterministically from the original input via the io.load_rna_ingest SSOT instead of
failing or forcing an S0 re-run. These tests pin the SSOT contract:
  - filtered input  -> pure load + a `counts` layer (no filtering);
  - raw input       -> deterministic barcode-rank knee cell-calling, returning the
                       pre-cell-call matrix (for SoupX) + diagnostics.
S0 and S1a both call this function, so they can never diverge on what "ingested RNA" is.
"""
import unittest
from unittest import mock

import anndata as ad
import numpy as np
import scipy.sparse as sp

from executor import io


def _toy(n=30, g=6):
    X = sp.csr_matrix((np.arange(n * g, dtype=float).reshape(n, g) + 1.0))
    return ad.AnnData(X=X)


class LoadRnaIngestSSOTTests(unittest.TestCase):
    def test_filtered_is_pure_load_plus_counts_layer(self):
        toy = _toy()
        with mock.patch.object(io, "load_rna", return_value=toy.copy()), \
             mock.patch.object(io, "detect_rna_format", return_value="10x_h5"):
            rna, raw_full, diag = io.load_rna_ingest("ignored", filtered_status="filtered")
        self.assertEqual(rna.n_obs, toy.n_obs)             # no filtering for filtered input
        self.assertEqual(rna.n_vars, toy.n_vars)
        self.assertIn("counts", rna.layers)
        self.assertEqual((rna.layers["counts"] != rna.X).nnz, 0)  # counts == X
        self.assertIsNone(raw_full)
        self.assertIsNone(diag)

    def test_raw_applies_deterministic_cell_calling(self):
        toy = _toy(n=60)
        with mock.patch.object(io, "load_rna", side_effect=lambda *a, **k: toy.copy()), \
             mock.patch.object(io, "detect_rna_format", return_value="10x_h5"):
            rna1, raw1, diag1 = io.load_rna_ingest("ignored", filtered_status="raw")
            rna2, raw2, diag2 = io.load_rna_ingest("ignored", filtered_status="raw")
        # Deterministic — identical cell set on repeat (so reconstruction matches S0).
        self.assertEqual(rna1.n_obs, rna2.n_obs)
        self.assertIn("counts", rna1.layers)
        self.assertIsNotNone(raw1)                         # pre-cell-call matrix (SoupX)
        self.assertEqual(raw1.n_obs, toy.n_obs)
        self.assertIsNotNone(diag1)                        # cell-calling diagnostics
        self.assertLessEqual(rna1.n_obs, toy.n_obs)        # cell-calling only ever drops


if __name__ == "__main__":
    unittest.main()
