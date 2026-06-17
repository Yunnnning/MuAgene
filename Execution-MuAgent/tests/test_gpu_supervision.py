"""GPU-utilisation liveness in the supervision daemon (classify rule + sstat parse)."""
import unittest
from unittest import mock

from execution_muagent import monitor


class ClassifyGpuTests(unittest.TestCase):
    def _ev(self, **kw):
        # Flat CPU + responsive fs + RUNNING would normally be confirmed_dead (rule 6).
        base = {"scheduler_state": "RUNNING", "error_markers": [], "child_storage_hang_ids": [],
                "filesystem_responsive": True, "cpu_delta": 0.0, "gpu_util": None}
        base.update(kw)
        return base

    def test_active_gpu_recovers_even_with_flat_cpu(self):
        verdict, reason = monitor.classify_investigation(self._ev(gpu_util=80.0))
        self.assertEqual(verdict, "recovered")
        self.assertIn("gpu_active", reason)

    def test_idle_gpu_does_not_save_a_dead_job(self):
        verdict, _ = monitor.classify_investigation(self._ev(gpu_util=0.0))
        self.assertEqual(verdict, "confirmed_dead")

    def test_absent_gpu_metric_falls_through_to_cpu_logic(self):
        verdict, _ = monitor.classify_investigation(self._ev(gpu_util=None))
        self.assertEqual(verdict, "confirmed_dead")


class ProbeGpuParseTests(unittest.TestCase):
    def test_parses_gpuutil_and_gpumem(self):
        with mock.patch.object(monitor.shutil, "which", return_value="/usr/bin/sstat"), \
             mock.patch.object(monitor, "_run_cmd",
                               return_value=(0, "cpu=00:01,gres/gpuutil=73,gres/gpumem=2048M", "")):
            res = monitor._probe_gpu("slurm", ["123"], 5)
        self.assertEqual(res["gpu_util"], 73.0)
        self.assertEqual(res["gpu_mem_mb"], 2048.0)

    def test_non_slurm_returns_none(self):
        self.assertEqual(monitor._probe_gpu("pbs", ["1"], 5), {"gpu_util": None, "gpu_mem_mb": None})

    def test_no_jobids_returns_none(self):
        self.assertEqual(monitor._probe_gpu("slurm", [], 5), {"gpu_util": None, "gpu_mem_mb": None})


if __name__ == "__main__":
    unittest.main()
