"""executor.compute device dispatch — cpu default, loud GPU failure, opt-in fallback."""
import os
import unittest
from unittest import mock

from executor import compute


class ComputeDeviceTests(unittest.TestCase):
    def setUp(self):
        self._env = dict(os.environ)

    def tearDown(self):
        os.environ.clear()
        os.environ.update(self._env)

    def test_default_is_cpu(self):
        os.environ.pop("PMA_DEVICE", None)
        self.assertEqual(compute.requested_device(), "cpu")
        self.assertFalse(compute.use_gpu())

    def test_gpu_requested_but_unavailable_raises_loud(self):
        os.environ["PMA_DEVICE"] = "gpu"
        os.environ.pop("PMA_DEVICE_FALLBACK", None)
        with mock.patch.object(compute, "gpu_capable", return_value=(False, "no rapids")):
            with self.assertRaises(compute.GpuUnavailableError):
                compute.use_gpu()

    def test_gpu_fallback_is_opt_in(self):
        os.environ["PMA_DEVICE"] = "gpu"
        os.environ["PMA_DEVICE_FALLBACK"] = "1"
        with mock.patch.object(compute, "gpu_capable", return_value=(False, "no rapids")):
            self.assertFalse(compute.use_gpu())  # warns + CPU, no raise

    def test_gpu_used_when_capable(self):
        os.environ["PMA_DEVICE"] = "gpu"
        with mock.patch.object(compute, "gpu_capable", return_value=(True, "ok")):
            self.assertTrue(compute.use_gpu())

    def test_device_used_labels(self):
        self.assertEqual(compute.device_used(True), "gpu")
        self.assertEqual(compute.device_used(False), "cpu")


if __name__ == "__main__":
    unittest.main()
