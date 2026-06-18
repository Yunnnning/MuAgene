"""Tests for s_handoff — per-sample post-QC handoff bundle (h5mu + manifest)."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import yaml

try:
    import anndata as ad
    import mudata as mu  # noqa: F401
    import numpy as np
    import scipy.sparse as sp
    _HAVE_DEPS = True
except Exception:
    _HAVE_DEPS = False

from executor.run_paths import RunPaths


def _params(branch: str, add_chr_prefix: bool, genome: str | None) -> dict:
    # provenance.get_value reads flat dotted keys wrapped as {"value": ...}.
    p = {
        "plan.workflow_branch": {"value": branch},
        "s2_atac_qc.add_chr_prefix": {"value": add_chr_prefix},
    }
    if genome is not None:
        p["ingest.genome_assembly"] = {"value": genome}
    return p


def _init_run(tmp: str, *, branch: str = "paired", add_chr_prefix: bool = True,
              genome: str | None = "mm10") -> RunPaths:
    paths = RunPaths(tmp)
    paths.ensure()
    paths.parameters_yaml.write_text(yaml.safe_dump(_params(branch, add_chr_prefix, genome)))
    return paths


def _write_rna(path: Path, n_obs: int = 6, n_var: int = 10) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    X = sp.csr_matrix(np.arange(n_obs * n_var, dtype="float32").reshape(n_obs, n_var))
    a = ad.AnnData(X=X)
    a.layers["counts"] = X.copy()
    a.obs_names = [f"cell{i}" for i in range(n_obs)]
    a.var_names = [f"gene{i}" for i in range(n_var)]
    a.write_h5ad(str(path))


def _write_atac(path: Path, n_obs: int = 6, n_var: int = 8) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    X = sp.csr_matrix(np.ones((n_obs, n_var), dtype="float32"))
    a = ad.AnnData(X=X)
    a.obs_names = [f"cell{i}" for i in range(n_obs)]
    a.var_names = [f"peak{i}" for i in range(n_var)]
    a.write_h5ad(str(path))


def _write_empty(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(X=sp.csr_matrix((0, 0))).write_h5ad(str(path))


def _run(paths: RunPaths, branch: str):
    from executor.stages import s_handoff
    return s_handoff.run(paths.run_dir, {}, workflow_branch=branch)


@unittest.skipUnless(_HAVE_DEPS, "anndata/mudata/numpy/scipy unavailable")
class SHandoffTests(unittest.TestCase):
    def test_paired_writes_h5mu_and_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="paired")
            _write_rna(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_atac(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))
            peaks = paths.artifact("s2_atac_qc", "peaks_macs3.bed")
            peaks.parent.mkdir(parents=True, exist_ok=True)
            peaks.write_text("chr1\t0\t100\n")
            paths.artifact("s2_atac_qc", "atac_fragments_cbf_chrnorm.tsv.gz").write_bytes(b"\x1f\x8b")

            result = _run(paths, "paired")

            h5mu = paths.deliv_results / f"post_qc_{paths.run_dir.name}.h5mu"
            self.assertTrue(h5mu.exists())
            mdata = mu.read_h5mu(str(h5mu))
            self.assertEqual(set(mdata.mod), {"rna", "atac"})
            self.assertEqual(mdata.mod["rna"].n_obs, 6)
            self.assertIn("counts", mdata.mod["rna"].layers)

            man = json.loads((paths.deliv_results / "post_qc_manifest.json").read_text())
            self.assertEqual(man["schema"], "muagene.post_qc_handoff/1")
            self.assertEqual(man["modality_branch"], "paired")
            self.assertEqual(man["genome_assembly"], "mm10")
            self.assertEqual(man["n_cells"], {"rna": 6, "atac": 6, "joint": 6})
            self.assertEqual(man["atac"]["peaks_bed"],
                             "internal/artifacts/s2_atac_qc/peaks_macs3.bed")
            self.assertEqual(man["atac"]["fragments_prepared"],
                             "internal/artifacts/s2_atac_qc/atac_fragments_cbf_chrnorm.tsv.gz")
            self.assertTrue(man["atac"]["add_chr_prefix"])
            self.assertEqual(man["atac"]["frag_chrom_convention"], "ucsc")
            self.assertTrue(result["post_qc_h5mu"].endswith(".h5mu"))

    def test_rna_only_drops_empty_atac_mod(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="rna_only", genome=None)
            _write_rna(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_empty(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))

            _run(paths, "rna_only")

            h5mu = paths.deliv_results / f"post_qc_{paths.run_dir.name}.h5mu"
            mdata = mu.read_h5mu(str(h5mu))
            self.assertEqual(set(mdata.mod), {"rna"})
            man = json.loads((paths.deliv_results / "post_qc_manifest.json").read_text())
            self.assertEqual(man["modality_branch"], "rna_only")
            self.assertIsNone(man["n_cells"]["atac"])
            self.assertIsNone(man["atac"]["peaks_bed"])
            self.assertIsNone(man["atac"]["frag_chrom_convention"])

    def test_atac_only_drops_empty_rna_mod(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="atac_only")
            _write_empty(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_atac(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))

            _run(paths, "atac_only")

            mdata = mu.read_h5mu(str(paths.deliv_results / f"post_qc_{paths.run_dir.name}.h5mu"))
            self.assertEqual(set(mdata.mod), {"atac"})
            man = json.loads((paths.deliv_results / "post_qc_manifest.json").read_text())
            self.assertEqual(man["modality_branch"], "atac_only")
            self.assertIsNone(man["n_cells"]["rna"])

    def test_raises_when_no_modality(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="paired")
            _write_empty(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_empty(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))
            with self.assertRaises(RuntimeError):
                _run(paths, "paired")


if __name__ == "__main__":
    unittest.main()
