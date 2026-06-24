"""Tests for _cleanup_process_intermediates — finish-phase S4–S8 working-file removal."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import yaml

from executor.cli import _cleanup_process_intermediates
from executor.run_paths import RunPaths
from executor.stage_progress import execute_done


def _init_run(tmp: str, *, branch: str = "paired") -> RunPaths:
    paths = RunPaths(tmp)
    paths.ensure()
    paths.parameters_yaml.write_text(
        yaml.safe_dump({"plan": {"workflow_branch": branch}})
    )
    return paths


def _touch(path: Path, content: bytes = b"placeholder") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


# The large RNA working h5ads and the ATAC sidecars/scratch finish-cleanup deletes.
_RNA_TARGETS = [
    ("s4_rna_norm", "rna_norm.h5ad"),
    ("s6_neighbors", "rna_neighbors.h5ad"),
    ("s7_clustering", "rna_clustered.h5ad"),
]
_ATAC_TARGETS = [
    ("s5_atac_spectral", "atac_spectral.h5ad"),
    ("s5_atac_spectral", "feature_matrix.npz"),
    ("s5_atac_spectral", "feature_names.tsv"),
    ("s5_atac_spectral", "feature_kind.txt"),
    ("s5_atac_spectral", "peak_matrix_s2peaks.h5ad"),
    ("s5_atac_spectral", "peak_matrix_user.h5ad"),
    ("s5_atac_spectral", "_s2_peaks_prepared.bed"),
    ("s5_atac_spectral", "_user_peaks_prepared.bed"),
    ("s7_clustering", "atac_leiden_labels.parquet"),
]
# Durable markers that MUST survive cleanup (status + DAG edges key off them).
_MARKERS = [
    ("s4_rna_norm", "norm_summary.json"),
    ("s5_atac_spectral", "spectral_summary.json"),
    ("s6_neighbors", "neighbors_summary.json"),
    ("s7_clustering", "clustering_summary.json"),
    ("s8_umap", "s8_done.txt"),
]


class CleanupProcessIntermediatesTests(unittest.TestCase):
    def test_deletes_s4_s7_working_files_paired(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="paired")
            targets = [_touch(paths.artifact(s, n)) for s, n in _RNA_TARGETS + _ATAC_TARGETS]
            markers = [_touch(paths.artifact(s, n), b"{}") for s, n in _MARKERS]
            deliverables = [
                _touch(paths.processed_h5mu),
                _touch(paths.run_manifest_json, b"{}"),
            ]

            deleted = _cleanup_process_intermediates(Path(tmp))

            for p in targets:
                self.assertFalse(p.exists(), f"Expected {p} to be deleted")
            self.assertEqual(len(deleted), len(targets))
            for p in markers + deliverables:
                self.assertTrue(p.exists(), f"Expected {p} to be preserved")

    def test_branch_aware_rna_only(self):
        """rna_only never wrote ATAC sidecars — only the RNA working files are removed."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="rna_only")
            rna = [_touch(paths.artifact(s, n)) for s, n in _RNA_TARGETS]
            # No S5/ATAC files created.

            deleted = _cleanup_process_intermediates(Path(tmp))

            for p in rna:
                self.assertFalse(p.exists())
            self.assertEqual(len(deleted), len(_RNA_TARGETS))
            self.assertFalse(any("s5_atac_spectral" in d for d in deleted),
                             "no S5 paths should be touched on rna_only")

    def test_branch_aware_atac_only(self):
        """atac_only: ATAC sidecars deleted, and the empty RNA stubs are removed too."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="atac_only")
            atac = [_touch(paths.artifact(s, n)) for s, n in _ATAC_TARGETS]
            rna_stubs = [_touch(paths.artifact(s, n)) for s, n in _RNA_TARGETS]

            deleted = _cleanup_process_intermediates(Path(tmp))

            for p in atac + rna_stubs:
                self.assertFalse(p.exists(), f"Expected {p} to be deleted")
            self.assertEqual(len(deleted), len(_ATAC_TARGETS) + len(_RNA_TARGETS))

    def test_preserves_deliverables_and_post_qc(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="paired")
            _touch(paths.artifact("s4_rna_norm", "rna_norm.h5ad"))
            kept = [
                _touch(paths.processed_h5mu),
                _touch(paths.run_manifest_json, b"{}"),
                _touch(paths.post_qc_h5mu),
            ]

            _cleanup_process_intermediates(Path(tmp))

            for p in kept:
                self.assertTrue(p.exists(), f"Expected {p} to be preserved")

    def test_returns_only_paths_that_existed(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="paired")
            _touch(paths.artifact("s4_rna_norm", "rna_norm.h5ad"))
            _touch(paths.artifact("s6_neighbors", "rna_neighbors.h5ad"))
            # The rest are absent.

            deleted = _cleanup_process_intermediates(Path(tmp))
            self.assertEqual(len(deleted), 2)

    def test_no_targets_present_returns_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            _init_run(tmp, branch="paired")
            self.assertEqual(_cleanup_process_intermediates(Path(tmp)), [])

    def test_stage_progress_still_done_after_cleanup(self):
        """Critical regression: deleting the working h5ads must NOT flip S4–S8 off.
        execute_done keys off the durable *_summary.json / s8_done.txt markers."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, branch="paired")
            # Markers + the deletable working files both present.
            for s, n in _MARKERS:
                _touch(paths.artifact(s, n), b"{}")
            for s, n in _RNA_TARGETS + _ATAC_TARGETS:
                _touch(paths.artifact(s, n))

            _cleanup_process_intermediates(Path(tmp))

            for stage in ("s4_rna_norm", "s5_atac_spectral", "s6_neighbors",
                          "s7_clustering", "s8_umap"):
                self.assertTrue(execute_done(paths, stage),
                                f"{stage} must still read as done after cleanup")


class ProcessMatrixUntrackedTests(unittest.TestCase):
    """The S4–S8 working files finish-cleanup deletes must NOT be declared Snakemake
    outputs, and the S6/S7/S8 input edges must depend on the durable *_summary.json
    markers — the structural fix that lets cleanup delete them without making any
    rule's declared output "missing" (which would re-run S4–S8 on the next submit).
    """

    _RULES = Path(__file__).resolve().parents[1] / "workflow" / "rules"

    @staticmethod
    def _code(text: str) -> str:
        """Drop comment lines so assertions test declarations, not prose."""
        return "\n".join(l for l in text.splitlines() if not l.strip().startswith("#"))

    def test_cleanup_targets_are_not_declared_outputs(self):
        names = [
            "rna_norm.h5ad", "rna_neighbors.h5ad", "rna_clustered.h5ad",
            "atac_spectral.h5ad", "feature_matrix.npz", "feature_names.tsv",
            "feature_kind.txt", "peak_matrix_s2peaks.h5ad", "peak_matrix_user.h5ad",
            "_s2_peaks_prepared.bed", "_user_peaks_prepared.bed",
            "atac_leiden_labels.parquet",
        ]
        joined = "\n".join(self._code(p.read_text())
                           for p in self._RULES.glob("s[4-8]*.smk"))
        for name in names:
            self.assertNotIn(name, joined,
                             f"{name} must not appear in any s4–s8 rule (declared output/edge)")

    def test_s6_s7_inputs_use_durable_markers_not_deletable_h5ads(self):
        s6 = self._code((self._RULES / "s6_neighbors.smk").read_text())
        self.assertIn("norm_summary.json", s6)
        self.assertNotIn("rna_norm.h5ad", s6)

        s7 = self._code((self._RULES / "s7_clustering.smk").read_text())
        self.assertIn("neighbors_summary.json", s7)
        self.assertNotIn("rna_neighbors.h5ad", s7)

        s8 = self._code((self._RULES / "s8_umap.smk").read_text())
        self.assertIn("clustering_summary.json", s8)
        self.assertNotIn("rna_clustered.h5ad", s8)


if __name__ == "__main__":
    unittest.main()
