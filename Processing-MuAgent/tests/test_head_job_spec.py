"""head_job.yaml carries `gpu_stages_present` = (_GPU_CAPABLE ∩ this run's stages) — the
flag Execution-MuAgent's execute-spec uses to gate the GPU env preflight. With the single
source `_GPU_CAPABLE` empty (preprocessing is CPU-only today) it is always False, so a
device=gpu run never pulls the GPU container for a stage that does not exist.
"""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import yaml

from executor import specs


def _flag(run_dir, target):
    out = specs.write_head_job_spec(run_dir, target)
    return yaml.safe_load(Path(out).read_text())["gpu_stages_present"]


class HeadJobGpuStagesPresentTests(unittest.TestCase):
    def test_false_when_gpu_capable_empty(self):
        # Real resources.smk: _GPU_CAPABLE is empty today.
        with tempfile.TemporaryDirectory() as tmp:
            self.assertFalse(_flag(tmp, "s3_doublets_execute"))

    def test_true_when_target_stage_is_gpu_capable(self):
        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch.object(specs, "_load_resources_smk",
                                   return_value=SimpleNamespace(_GPU_CAPABLE={"s3_doublets"})):
                self.assertTrue(_flag(tmp, "s3_doublets_execute"))

    def test_false_when_a_different_stage_is_gpu_capable(self):
        # The flag tracks the *target* stage, not merely "some GPU stage exists".
        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch.object(specs, "_load_resources_smk",
                                   return_value=SimpleNamespace(_GPU_CAPABLE={"s6_neighbors"})):
                self.assertFalse(_flag(tmp, "s3_doublets_execute"))


if __name__ == "__main__":
    unittest.main()
