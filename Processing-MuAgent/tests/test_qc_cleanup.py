"""Tests for _cleanup_qc_intermediates — post-approval h5ad removal."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import yaml

from executor.cli import _cleanup_qc_intermediates
from executor.run_paths import RunPaths


def _init_run(tmp: str) -> RunPaths:
    paths = RunPaths(tmp)
    paths.ensure()
    paths.parameters_yaml.write_text(
        yaml.safe_dump({"plan": {"workflow_branch": "paired"}})
    )
    return paths


def _touch(path: Path, content: bytes = b"placeholder") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


class CleanupQCIntermediatesTests(unittest.TestCase):
    def test_deletes_target_h5ads(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            rna_qc   = _touch(paths.artifact("s1_rna_qc",  "rna_qc.h5ad"))
            atac_qc  = _touch(paths.artifact("s2_atac_qc", "atac_qc.h5ad"))
            atac_snap = _touch(paths.artifact("s2_atac_qc", "atac_snap.h5ad"))
            atac_snap_explore = _touch(paths.artifact("qc_explore", "atac_snap_explore.h5ad"))

            deleted = _cleanup_qc_intermediates(Path(tmp))

            self.assertFalse(rna_qc.exists())
            self.assertFalse(atac_qc.exists())
            self.assertFalse(atac_snap.exists())
            self.assertFalse(atac_snap_explore.exists())
            self.assertEqual(len(deleted), 4)

    def test_returns_only_paths_that_existed(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            # Only create two of the three targets
            _touch(paths.artifact("s1_rna_qc",  "rna_qc.h5ad"))
            _touch(paths.artifact("s2_atac_qc", "atac_qc.h5ad"))
            # atac_snap.h5ad is absent (already cleaned or never created)

            deleted = _cleanup_qc_intermediates(Path(tmp))
            self.assertEqual(len(deleted), 2)

    def test_preserves_qc_summary_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            _touch(paths.artifact("s1_rna_qc",  "rna_qc.h5ad"))
            _touch(paths.artifact("s2_atac_qc", "atac_qc.h5ad"))
            _touch(paths.artifact("s2_atac_qc", "atac_snap.h5ad"))
            s1_json  = _touch(paths.artifact("s1_rna_qc",  "qc_summary.json"), b"{}")
            s2_json  = _touch(paths.artifact("s2_atac_qc", "qc_summary.json"), b"{}")

            _cleanup_qc_intermediates(Path(tmp))

            self.assertTrue(s1_json.exists())
            self.assertTrue(s2_json.exists())

    def test_preserves_parquets_and_cbf_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            _touch(paths.artifact("s1_rna_qc", "rna_qc.h5ad"))
            kept = [
                _touch(paths.artifact("s1_rna_qc", "qc_metrics_pre.parquet"),  b"PAR1\x00PAR1"),
                _touch(paths.artifact("s1_rna_qc", "qc_metrics_post.parquet"), b"PAR1\x00PAR1"),
                _touch(paths.artifact("s2_atac_qc", "atac_fragments_cbf.tsv.gz"), b"data"),
            ]

            _cleanup_qc_intermediates(Path(tmp))

            for p in kept:
                self.assertTrue(p.exists(), f"Expected {p} to be preserved")

    def test_preserves_s1a_and_s3_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            kept = [
                _touch(paths.artifact("s1a_ambient", "rna_decontaminated.h5ad")),
                _touch(paths.artifact("s3_doublets", "rna_post_doublet.h5ad")),
                _touch(paths.artifact("s3_doublets", "atac_post_doublet.h5ad")),
                _touch(paths.artifact("s3_doublets", "calls.parquet"), b"PAR1\x00PAR1"),
            ]

            _cleanup_qc_intermediates(Path(tmp))

            for p in kept:
                self.assertTrue(p.exists(), f"Expected {p} to be preserved")

    def test_preserves_qc_explore_metric_parquets(self):
        """The per-cell QC metric parquets (under qc_explore/) must survive cleanup
        so a post-approval `revise` can re-derive thresholds without a heavy reload."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            _touch(paths.artifact("s1_rna_qc", "rna_qc.h5ad"))
            _touch(paths.artifact("s2_atac_qc", "atac_qc.h5ad"))
            kept = [
                _touch(paths.artifact("qc_explore", "rna_qc_metrics.parquet"), b"PAR1\x00PAR1"),
                _touch(paths.artifact("qc_explore", "atac_qc_metrics.parquet"), b"PAR1\x00PAR1"),
                _touch(paths.artifact("qc_explore", "qc_explore.json"), b"{}"),
            ]

            _cleanup_qc_intermediates(Path(tmp))

            for p in kept:
                self.assertTrue(p.exists(), f"Expected {p} to be preserved")

    def test_no_targets_present_returns_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            _init_run(tmp)
            deleted = _cleanup_qc_intermediates(Path(tmp))
            self.assertEqual(deleted, [])


if __name__ == "__main__":
    unittest.main()
