"""Tests for central deliverables/figures/ layout and lazy scaffold."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from executor.run_paths import RunPaths


class DeliverablesFigureLayoutTests(unittest.TestCase):
    def test_ensure_creates_plan_only_not_later_phases(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = RunPaths(tmp)
            paths.ensure()
            self.assertTrue(paths.deliv_plan.is_dir())
            self.assertFalse(paths.deliv_figures.exists())
            self.assertFalse(paths.deliv_checkpoints.exists())
            self.assertFalse(paths.deliv_results.exists())

    def test_deliv_figures_path_uses_central_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = RunPaths(tmp)
            paths.ensure()
            expected = paths.deliv_figures / "s1_rna_qc_violin_pre.png"
            self.assertEqual(paths.deliv_figures_path("s1_rna_qc_violin_pre"), expected)

    def test_migrate_moves_legacy_figures_into_central_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = RunPaths(tmp)
            paths.ensure()
            legacy_qc = paths.deliverables / "checkpoint" / "qc_review" / "figures"
            legacy_qc.mkdir(parents=True, exist_ok=True)
            flat_png = legacy_qc / "post_qc_review_cell_counts.png"
            flat_png.write_bytes(b"PNG")
            paths.deliv_checkpoints.mkdir(parents=True, exist_ok=True)
            paths.deliv_qc_review.mkdir(parents=True, exist_ok=True)
            paths.qc_review_summary_md.write_text("# QC review\n")

            moved = paths.migrate_legacy_figures_to_central()

            self.assertFalse(flat_png.exists())
            self.assertTrue(paths.deliv_figures_path("post_qc_review_cell_counts").is_file())
            self.assertEqual(len(moved), 1)
            self.assertTrue(paths.qc_review_summary_md.is_file())

    def test_resolve_figure_falls_back_to_legacy_qc_review_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = RunPaths(tmp)
            paths.ensure()
            legacy = paths.deliverables / "figure" / "s2_atac_qc_frip_histogram.png"
            legacy.parent.mkdir(parents=True, exist_ok=True)
            legacy.write_bytes(b"PNG")

            resolved = paths.resolve_figure("s2_atac_qc_frip_histogram")

            self.assertEqual(resolved, legacy)

    def test_resolve_figure_falls_back_to_legacy_post_run_umap(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = RunPaths(tmp)
            paths.ensure()
            legacy = paths.deliverables / "post_run" / "s8_umap_rna_by_leiden.png"
            legacy.parent.mkdir(parents=True, exist_ok=True)
            legacy.write_bytes(b"PNG")

            resolved = paths.resolve_figure("s8_umap_rna_by_leiden")

            self.assertEqual(resolved, legacy)


if __name__ == "__main__":
    unittest.main()
