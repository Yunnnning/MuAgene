"""Regression test: revising a QC stage invalidates the right artifacts + gate.

Re-running a QC stage after `revise` only works if the stale downstream
artifacts AND the post_qc_review_propose gate outputs are deleted. Otherwise
Snakemake reports "Nothing to be done" (the terminal target looks satisfied) and
silently skips the re-run. `_invalidate_qc_downstream` must delete exactly the
right set, preserve the expensive chr-norm cache, and be a no-op for non-QC
stages / missing files.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from executor.cli import _invalidate_qc_downstream
from executor.run_paths import RunPaths


class QcInvalidationTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.run_dir = Path(self._tmp.name)
        self.paths = RunPaths(self.run_dir)
        self.paths.ensure()

    def tearDown(self):
        self._tmp.cleanup()

    def _touch(self, p: Path):
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("x")
        return p

    def _make_full_qc_tree(self):
        rp = self.paths
        files = {
            "s1_rna_qc": [rp.artifact("s1_rna_qc", "rna_qc.h5ad"),
                          rp.artifact("s1_rna_qc", "qc_summary.json")],
            "s2_atac_qc": [rp.artifact("s2_atac_qc", "atac_qc.h5ad"),
                           rp.artifact("s2_atac_qc", "qc_summary.json")],
            "cache": [rp.artifact("s2_atac_qc", "atac_fragments_cbf_chrnorm.tsv.gz")],
            "s3": [rp.artifact("s3_doublets", "rna_post_doublet.h5ad"),
                   rp.artifact("s3_doublets", "atac_post_doublet.h5ad"),
                   rp.artifact("s3_doublets", "calls.parquet"),
                   rp.artifact("s3_doublets", "joint_barcodes.txt"),
                   rp.artifact("s3_doublets", "overlap_summary.json")],
            "gate": [rp.proposal("post_qc_review"), rp.awaiting_sentinel("post_qc_review"),
                     rp.qc_review_summary_md, rp.qc_summary_html],
        }
        for group in files.values():
            for p in group:
                self._touch(p)
        return files

    def test_s2_revision_clears_s2_s3_and_gate_but_keeps_cache_and_s1(self):
        f = self._make_full_qc_tree()
        deleted = _invalidate_qc_downstream(self.run_dir, "s2_atac_qc")
        # S2 + S3 + gate gone
        for p in f["s2_atac_qc"] + f["s3"] + f["gate"]:
            self.assertFalse(p.exists(), f"should be deleted: {p}")
        # Expensive chr-norm cache preserved
        for p in f["cache"]:
            self.assertTrue(p.exists(), f"cache must be preserved: {p}")
        # S1 (not downstream of S2) untouched
        for p in f["s1_rna_qc"]:
            self.assertTrue(p.exists(), f"S1 must be untouched on S2 revision: {p}")
        self.assertTrue(any("post_qc_review" in d for d in deleted))

    def test_s3_revision_clears_only_s3_and_gate(self):
        f = self._make_full_qc_tree()
        _invalidate_qc_downstream(self.run_dir, "s3_doublets")
        for p in f["s3"] + f["gate"]:
            self.assertFalse(p.exists(), f"should be deleted: {p}")
        # S2 NOT downstream of S3 — untouched
        for p in f["s2_atac_qc"]:
            self.assertTrue(p.exists(), f"S2 must be untouched on S3 revision: {p}")

    def test_non_qc_stage_is_noop(self):
        self._make_full_qc_tree()
        self.assertEqual(_invalidate_qc_downstream(self.run_dir, "plan_review"), [])

    def test_missing_files_noop_no_error(self):
        # Nothing created — must not raise, returns [].
        self.assertEqual(_invalidate_qc_downstream(self.run_dir, "s2_atac_qc"), [])


if __name__ == "__main__":
    unittest.main()
