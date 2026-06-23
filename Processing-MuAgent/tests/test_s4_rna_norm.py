"""S4 reads post-QC RNA from the canonical handoff h5mu, with a legacy s3-h5ad fallback."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import yaml

try:
    import anndata as ad
    import mudata as mu
    import numpy as np
    import scipy.sparse as sp
    _HAVE = True
except Exception:
    _HAVE = False

from executor.run_paths import RunPaths

_PLAN = {"stages": {"s4_rna_norm": {"parameters": {
    "target_sum": {"value": 1e4},
    "hvg_flavor": {"value": "seurat_v3"},
    "hvg_n_top_genes": {"value": 2000},
}}}}


def _rna(n_obs: int = 50, n_var: int = 300):
    rng = np.random.default_rng(0)
    X = sp.csr_matrix(rng.integers(0, 6, size=(n_obs, n_var)).astype("float32"))
    a = ad.AnnData(X=X.copy())
    a.layers["counts"] = X.copy()
    a.obs_names = [f"cell{i}" for i in range(n_obs)]
    a.var_names = [f"gene{i}" for i in range(n_var)]
    return a


def _init(tmp: str, branch: str = "rna_only") -> RunPaths:
    paths = RunPaths(tmp)
    paths.ensure()
    paths.parameters_yaml.write_text(yaml.safe_dump({"plan.workflow_branch": {"value": branch}}))
    return paths


@unittest.skipUnless(_HAVE, "anndata/mudata/numpy/scipy unavailable")
class S4ReadsH5muTests(unittest.TestCase):
    def test_reads_rna_from_h5mu_without_s3_h5ad(self):
        from executor.stages import s4_rna_norm
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init(tmp, "rna_only")
            paths.deliv_qc.mkdir(parents=True, exist_ok=True)
            mu.MuData({"rna": _rna()}).write(str(paths.post_qc_h5mu), compression="gzip")
            # The transient s3 h5ad is absent — S4 must work from the h5mu alone.
            self.assertFalse(paths.artifact("s3_doublets", "rna_post_doublet.h5ad").exists())

            s4_rna_norm.run(paths.run_dir, _PLAN)

            out = paths.artifact("s4_rna_norm", "rna_norm.h5ad")
            self.assertTrue(out.exists())
            self.assertIn("highly_variable", ad.read_h5ad(out).var)

    def test_legacy_fallback_to_s3_h5ad(self):
        from executor.stages import s4_rna_norm
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init(tmp, "rna_only")
            legacy = paths.artifact("s3_doublets", "rna_post_doublet.h5ad")
            legacy.parent.mkdir(parents=True, exist_ok=True)
            _rna().write_h5ad(str(legacy))
            self.assertFalse(paths.post_qc_h5mu.exists())  # no h5mu → fall back

            s4_rna_norm.run(paths.run_dir, _PLAN)
            self.assertTrue(paths.artifact("s4_rna_norm", "rna_norm.h5ad").exists())

    def test_raises_when_neither_source_present(self):
        from executor.stages import s4_rna_norm
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init(tmp, "rna_only")
            with self.assertRaises(FileNotFoundError):
                s4_rna_norm.run(paths.run_dir, _PLAN)


if __name__ == "__main__":
    unittest.main()
