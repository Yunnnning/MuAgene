"""Tests for qc_handoff — per-sample post-QC handoff bundle (h5mu + manifest).

The ATAC fixture is written with SnapATAC2's native (Blosc-compressed) writer via
``snap.pp.import_fragments`` — the way S3 actually writes ``atac_post_doublet.h5ad``.
A plain-anndata fixture would NOT exercise ``qc_handoff._load_atac_mod``'s snap.read
path and so would not have caught the ATAC-dropped-from-h5mu regression this guards.
"""
from __future__ import annotations

import gzip
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
    import snapatac2 as snap
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


def _write_atac_snap(path: Path, n_obs: int = 6) -> None:
    """Write a SnapATAC2-native post-doublet ATAC h5ad (Blosc-compressed) the way S3
    does — via import_fragments — so qc_handoff reads it through snap.read."""
    path.parent.mkdir(parents=True, exist_ok=True)
    barcodes = [f"cell{i}" for i in range(n_obs)]
    rng = np.random.default_rng(0)
    rows: list[str] = []
    for bc in barcodes:
        for s in sorted(int(x) for x in rng.integers(0, 9000, size=50)):
            rows.append(f"chr1\t{s}\t{s + 50}\t{bc}\t1")
    rows.sort(key=lambda r: int(r.split("\t")[1]))
    frag = path.parent / "_frag_fixture.tsv.gz"
    with gzip.open(frag, "wt") as fh:
        fh.write("\n".join(rows) + "\n")
    data = snap.pp.import_fragments(
        str(frag), chrom_sizes={"chr1": 10_000}, is_paired=True,
        file=str(path), min_num_fragments=1, sorted_by_barcode=False,
    )
    try:
        data.close()
    except Exception:
        pass


def _write_empty(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(X=sp.csr_matrix((0, 0))).write_h5ad(str(path))


def _run(paths: RunPaths, branch: str):
    from executor.stages import qc_handoff
    return qc_handoff.run(paths.run_dir, {}, workflow_branch=branch)


@unittest.skipUnless(_HAVE_DEPS, "anndata/mudata/numpy/scipy/snapatac2 unavailable")
class QcHandoffTests(unittest.TestCase):
    def test_paired_writes_h5mu_with_both_mods(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="paired")
            _write_rna(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_atac_snap(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))
            peaks = paths.artifact("s2_atac_qc", "peaks_macs3.bed")
            peaks.parent.mkdir(parents=True, exist_ok=True)
            peaks.write_text("chr1\t0\t100\n")
            paths.artifact("s2_atac_qc", "atac_fragments_cbf_chrnorm.tsv.gz").write_bytes(b"\x1f\x8b")

            result = _run(paths, "paired")

            h5mu = paths.post_qc_h5mu
            self.assertTrue(h5mu.exists())
            mdata = mu.read_h5mu(str(h5mu))
            # Regression guard: ATAC must be present on a paired run (was silently
            # dropped because anndata.read_h5ad cannot decode the snap Blosc matrix).
            self.assertEqual(set(mdata.mod), {"rna", "atac"})
            self.assertEqual(mdata.mod["rna"].n_obs, 6)
            self.assertIn("counts", mdata.mod["rna"].layers)
            # Lean ATAC mod carries fragments + chrom sizes for a downstream re-tile.
            atac = mdata.mod["atac"]
            self.assertEqual(atac.n_obs, 6)
            self.assertIn("fragment_paired", atac.obsm)
            self.assertIn("reference_sequences", atac.uns)

            man = json.loads(paths.post_qc_manifest_json.read_text())
            self.assertEqual(man["schema"], "muagene.post_qc_handoff/1")
            self.assertEqual(man["modality_branch"], "paired")
            self.assertEqual(man["genome_assembly"], "mm10")
            # n_cells.atac must be a real count, not null.
            self.assertEqual(man["n_cells"]["rna"], 6)
            self.assertEqual(man["n_cells"]["atac"], 6)
            self.assertIsNotNone(man["n_cells"]["atac"])
            self.assertEqual(man["atac"]["peaks_bed"],
                             "internal/artifacts/s2_atac_qc/peaks_macs3.bed")
            self.assertEqual(man["atac"]["fragments_prepared"],
                             "internal/artifacts/s2_atac_qc/atac_fragments_cbf_chrnorm.tsv.gz")
            self.assertTrue(man["atac"]["add_chr_prefix"])
            self.assertEqual(man["atac"]["frag_chrom_convention"], "ucsc")
            self.assertIn("deliverables/qc/", result["post_qc_h5mu"])
            # Dedup: the redundant post-doublet h5ads are deleted after the bundle.
            self.assertFalse(paths.artifact("s3_doublets", "rna_post_doublet.h5ad").exists())
            self.assertFalse(paths.artifact("s3_doublets", "atac_post_doublet.h5ad").exists())
            self.assertEqual(len(result["deleted_post_doublet_h5ads"]), 2)

    def test_atac_mod_is_retileable(self):
        """The lean ATAC mod, extracted from the h5mu, must rebuild its tile matrix —
        this is what the deferred S5-from-h5mu refactor will rely on."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="atac_only")
            _write_empty(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_atac_snap(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))

            _run(paths, "atac_only")

            atac = mu.read_h5mu(str(paths.post_qc_h5mu)).mod["atac"]
            standalone = Path(tmp) / "atac_extracted.h5ad"
            atac.write_h5ad(str(standalone), compression="gzip")
            sa = snap.read(str(standalone), backed=None)  # in-memory: dims can grow
            snap.pp.add_tile_matrix(sa, bin_size=500)
            self.assertGreater(sa.n_vars, 0)  # tile matrix rebuilt from fragments

    def test_rna_only_drops_empty_atac_mod(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="rna_only", genome=None)
            _write_rna(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_empty(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))

            _run(paths, "rna_only")

            mdata = mu.read_h5mu(str(paths.post_qc_h5mu))
            self.assertEqual(set(mdata.mod), {"rna"})
            man = json.loads(paths.post_qc_manifest_json.read_text())
            self.assertEqual(man["modality_branch"], "rna_only")
            self.assertIsNone(man["n_cells"]["atac"])
            self.assertIsNone(man["atac"]["peaks_bed"])
            self.assertIsNone(man["atac"]["frag_chrom_convention"])

    def test_atac_only_drops_empty_rna_mod(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="atac_only")
            _write_empty(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_atac_snap(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))

            _run(paths, "atac_only")

            mdata = mu.read_h5mu(str(paths.post_qc_h5mu))
            self.assertEqual(set(mdata.mod), {"atac"})
            man = json.loads(paths.post_qc_manifest_json.read_text())
            self.assertEqual(man["modality_branch"], "atac_only")
            self.assertIsNone(man["n_cells"]["rna"])

    def test_raises_when_no_modality(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="paired")
            _write_empty(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_empty(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))
            with self.assertRaises(RuntimeError):
                _run(paths, "paired")

    def test_deletes_h5ads_but_keeps_s3_markers(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="paired")
            _write_rna(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_atac_snap(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))
            # S3 markers that must survive the dedup deletion.
            calls = paths.artifact("s3_doublets", "calls.parquet"); calls.write_bytes(b"PARQUET")
            jb = paths.artifact("s3_doublets", "joint_barcodes.txt"); jb.write_text("cell0\n")
            ov = paths.artifact("s3_doublets", "overlap_summary.json"); ov.write_text("{}")

            _run(paths, "paired")

            self.assertFalse(paths.artifact("s3_doublets", "rna_post_doublet.h5ad").exists())
            self.assertFalse(paths.artifact("s3_doublets", "atac_post_doublet.h5ad").exists())
            self.assertTrue(calls.exists())
            self.assertTrue(jb.exists())
            self.assertTrue(ov.exists())

    def test_atac_only_deletes_h5ads(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="atac_only")
            _write_empty(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_atac_snap(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))
            _run(paths, "atac_only")
            self.assertFalse(paths.artifact("s3_doublets", "atac_post_doublet.h5ad").exists())
            self.assertFalse(paths.artifact("s3_doublets", "rna_post_doublet.h5ad").exists())

    def test_s3_declares_only_durable_marker(self):
        """Durability tripwire: S3 declares calls.parquet (durable) and NOT the
        post-doublet h5ads, so qc_handoff deleting them never triggers an S3 re-run."""
        from executor import specs
        outs = " ".join(specs._STAGE_IO["s3_doublets"]["outputs"].values())
        self.assertIn("calls.parquet", outs)
        self.assertNotIn("post_doublet", outs)

    def test_raises_when_branch_expects_atac_but_missing(self):
        """Loud failure: a paired run with RNA but an empty ATAC placeholder must
        raise rather than silently emit an RNA-only bundle (the old bug)."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="paired")
            _write_rna(paths.artifact("s3_doublets", "rna_post_doublet.h5ad"))
            _write_empty(paths.artifact("s3_doublets", "atac_post_doublet.h5ad"))
            with self.assertRaises(RuntimeError):
                _run(paths, "paired")


if __name__ == "__main__":
    unittest.main()
