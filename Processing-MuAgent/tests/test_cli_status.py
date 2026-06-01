import tempfile
import unittest

from executor import cli
from executor.run_paths import RunPaths


class CliStatusTests(unittest.TestCase):
    def test_checkpoint_aliases_map_to_internal_stage_names(self):
        self.assertEqual(cli._canonical_stage("qc_review"), "post_qc_review")
        self.assertEqual(cli._canonical_stage("resolution_review"), "s7_clustering")
        self.assertEqual(cli._display_stage("post_qc_review"), "qc_review")
        self.assertEqual(cli._display_stage("s7_clustering"), "resolution_review")

    def test_stage_states_include_pipeline_substeps_and_resolution_review(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = RunPaths(tmp)
            paths.ensure()
            paths.awaiting_sentinel("post_qc_review").write_text("")

            rows = cli._stage_states(paths)
            labels = [r[0] for r in rows]
            tasks = [r[1] for r in rows]
            states = {r[0]: r[2] for r in rows}

            self.assertEqual(states["qc_review"], "awaiting_approval")
            self.assertIn("resolution_review", labels)
            self.assertIn("S1a", labels)
            self.assertIn("Ambient RNA correction", tasks)
            self.assertIn("Resolution review", tasks)


if __name__ == "__main__":
    unittest.main()
