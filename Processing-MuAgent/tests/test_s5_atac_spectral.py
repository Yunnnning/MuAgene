"""S5 builds its SnapATAC2 working object from the post-QC h5mu, with a legacy fallback.

Focuses on `_build_atac_working` (the dedup change): the rest of S5 (select_features /
spectral / peak export) is exercised end-to-end on real data in the whelanC57A run. The
ATAC fixture is snap-native (import_fragments), as S3 writes it.
"""
from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

import yaml

try:
    import anndata as ad
    import mudata as mu  # noqa: F401
    import numpy as np
    import scipy.sparse as sp
    import snapatac2 as snap
    _HAVE = True
except Exception:
    _HAVE = False

from executor.run_paths import RunPaths


def _write_atac_snap(path: Path, n_obs: int = 6) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    bcs = [f"cell{i}" for i in range(n_obs)]
    rng = np.random.default_rng(0)
    rows: list[str] = []
    for bc in bcs:
        for s in sorted(int(x) for x in rng.integers(0, 9000, size=50)):
            rows.append(f"chr1\t{s}\t{s + 50}\t{bc}\t1")
    rows.sort(key=lambda r: int(r.split("\t")[1]))
    frag = path.parent / "_frag.tsv.gz"
    with gzip.open(frag, "wt") as fh:
        fh.write("\n".join(rows) + "\n")
    d = snap.pp.import_fragments(
        str(frag), chrom_sizes={"chr1": 10_000}, is_paired=True,
        file=str(path), min_num_fragments=1, sorted_by_barcode=False,
    )
    try:
        d.close()
    except Exception:
        pass


def _init(tmp: str) -> RunPaths:
    paths = RunPaths(tmp)
    paths.ensure()
    paths.parameters_yaml.write_text(yaml.safe_dump({
        "plan.workflow_branch": {"value": "atac_only"},
        "s2_atac_qc.add_chr_prefix": {"value": False},
    }))
    return paths


@unittest.skipUnless(_HAVE, "anndata/mudata/numpy/scipy/snapatac2 unavailable")
class S5BuildAtacWorkingTests(unittest.TestCase):
    def test_builds_snap_native_working_from_h5mu(self):
        from executor.stages import qc_handoff, s5_atac_spectral as s5
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init(tmp)
            _write_atac_snap(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))
            ad.AnnData(X=sp.csr_matrix((0, 0))).write_h5ad(
                str(paths.artifact("s3_doublets", "rna_post_doublet.h5ad")))
            # Produce the post-QC h5mu (atac mod = lean fragments); deletes the s3 h5ads.
            qc_handoff.run(paths.run_dir, {}, workflow_branch="atac_only")
            self.assertTrue(paths.post_qc_h5mu.exists())
            self.assertFalse(paths.artifact("s3_doublets", "atac_post_doublet.h5ad").exists())

            art = Path(tmp) / "internal" / "artifacts" / "s5_atac_spectral"
            art.mkdir(parents=True, exist_ok=True)
            dst = art / "atac_spectral.h5ad"
            adata = s5._build_atac_working(paths.run_dir, dst)
            # Regression: barcodes must survive — snap.AnnData drops a pandas obs
            # index and would assign integer labels ('0','1',...), breaking S8's
            # RNA<->ATAC barcode intersection. Must match the cells (not integers).
            self.assertEqual(set(adata.obs_names), {f"cell{i}" for i in range(6)})
            self.assertFalse(all(str(b).isdigit() for b in adata.obs_names))
            # snap-native + re-tileable from the h5mu fragments (no s3 h5ad on disk).
            snap.pp.add_tile_matrix(adata, bin_size=500)
            self.assertEqual(adata.n_obs, 6)
            self.assertGreater(adata.n_vars, 0)
            try:
                adata.close()
            except Exception:
                pass

    def test_legacy_fallback_to_s3_h5ad(self):
        from executor.stages import s5_atac_spectral as s5
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init(tmp)
            _write_atac_snap(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))
            self.assertFalse(paths.post_qc_h5mu.exists())  # no h5mu → legacy copy path

            art = Path(tmp) / "internal" / "artifacts" / "s5_atac_spectral"
            art.mkdir(parents=True, exist_ok=True)
            dst = art / "atac_spectral.h5ad"
            adata = s5._build_atac_working(paths.run_dir, dst)
            self.assertEqual(adata.n_obs, 6)
            self.assertEqual(set(adata.obs_names), {f"cell{i}" for i in range(6)})
            try:
                adata.close()
            except Exception:
                pass


if __name__ == "__main__":
    unittest.main()
