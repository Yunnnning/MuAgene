"""write_stage_specs refreshes the current branch's per-stage specs AND prunes orphan
specs left behind by a renamed/branch-dropped stage (e.g. s_handoff -> qc_handoff), while
always preserving head_job.yaml. Guards the stale-spec-on-resubmit gap: a resubmit must not
leave Execution-MuAgent's monitor verifying a stage that no longer exists.
"""
import tempfile
import unittest
from pathlib import Path

from executor import specs
from executor.run_paths import RunPaths


class StageSpecPruneTests(unittest.TestCase):
    def test_prunes_orphan_and_preserves_head_job(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = RunPaths(Path(tmp))
            paths.stage_meta_dir.mkdir(parents=True, exist_ok=True)
            # An orphan left by a renamed stage, plus the separately-written head job.
            (paths.stage_meta_dir / "s_handoff.yaml").write_text("stage: s_handoff\n")
            (paths.stage_meta_dir / "head_job.yaml").write_text("stage: head_job\n")

            written = specs.write_stage_specs(tmp, "paired")
            names = {p.name for p in paths.stage_meta_dir.glob("*.yaml")}

            self.assertFalse(
                (paths.stage_meta_dir / "s_handoff.yaml").exists(),
                "orphan s_handoff.yaml should be pruned on respec",
            )
            self.assertIn("qc_handoff.yaml", names, "current handoff spec should be written")
            self.assertIn("head_job.yaml", names, "head_job.yaml must always be preserved")
            # Exactly the freshly-written branch specs + head_job.yaml remain.
            self.assertEqual(names, {p.name for p in written} | {"head_job.yaml"})

    def test_branch_change_prunes_dropped_stage(self):
        # rna_only has no ATAC stages; a prior paired-run spec must not survive.
        with tempfile.TemporaryDirectory() as tmp:
            paths = RunPaths(Path(tmp))
            paths.stage_meta_dir.mkdir(parents=True, exist_ok=True)
            (paths.stage_meta_dir / "s2_atac_qc.yaml").write_text("stage: s2_atac_qc\n")

            specs.write_stage_specs(tmp, "rna_only")
            names = {p.name for p in paths.stage_meta_dir.glob("*.yaml")}

            self.assertNotIn("s2_atac_qc.yaml", names,
                             "ATAC spec should be pruned for an rna_only branch")

    def test_qc_handoff_spec_peaks_bed_branch_aware(self):
        """peaks_bed is declared for ATAC branches only — matches qc_handoff.smk outputs."""
        import yaml
        with tempfile.TemporaryDirectory() as tmp:
            specs.write_stage_specs(tmp, "paired")
            paired = yaml.safe_load((RunPaths(Path(tmp)).stage_meta("qc_handoff")).read_text())
            self.assertIn("peaks_bed", paired["outputs"])

            specs.write_stage_specs(tmp, "rna_only")
            rna_only = yaml.safe_load((RunPaths(Path(tmp)).stage_meta("qc_handoff")).read_text())
            self.assertNotIn("peaks_bed", rna_only["outputs"])


if __name__ == "__main__":
    unittest.main()
