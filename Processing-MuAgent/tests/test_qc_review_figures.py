"""Tests for QC review figure subdirectory layout."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from executor.run_paths import RunPaths


class QCReviewFiguresLayoutTests(unittest.TestCase):
    def test_deliv_qc_figure_uses_figures_subdir(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = RunPaths(tmp)
            paths.ensure()
            expected = paths.deliv_qc_review / "figures" / "s1_rna_qc_violin_pre.png"
            self.assertEqual(paths.deliv_qc_figure("s1_rna_qc_violin_pre"), expected)

    def test_migrate_moves_flat_figures_into_subdir(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = RunPaths(tmp)
            paths.ensure()
            flat_png = paths.deliv_qc_review / "post_qc_review_cell_counts.png"
            flat_png.write_bytes(b"PNG")
            paths.qc_review_summary_md.write_text("# QC review\n")

            moved = paths.migrate_qc_figures_to_subdir()

            self.assertFalse(flat_png.exists())
            self.assertTrue(paths.deliv_qc_figure("post_qc_review_cell_counts").is_file())
            self.assertEqual(len(moved), 1)
            self.assertTrue(paths.qc_review_summary_md.is_file())

    def test_resolve_qc_figure_falls_back_to_legacy_flat_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = RunPaths(tmp)
            paths.ensure()
            legacy = paths.deliv_qc_review / "s2_atac_qc_frip_histogram.png"
            legacy.write_bytes(b"PNG")

            resolved = paths.resolve_qc_figure("s2_atac_qc_frip_histogram")

            self.assertEqual(resolved, legacy)


if __name__ == "__main__":
    unittest.main()
