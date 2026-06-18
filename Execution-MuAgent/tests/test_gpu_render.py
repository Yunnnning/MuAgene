"""render_submission_script + load_site_config carry the device/GPU routing vars.

The head job stays CPU; these exports only steer the GPU-capable child jobs that
Snakemake submits via the profile submit scripts.
"""
from pathlib import Path
import tempfile
import unittest

from execution_muagent.monitor import (
    SiteConfig,
    StageSpec,
    load_site_config,
    render_submission_script,
)


def _spec(stage="s3_doublets"):
    return StageSpec(
        schema_version="1",
        stage=stage,
        science_description="test",
        resources={"cpus": 2, "mem_mb": 8000, "walltime_min": 60},
        inputs={},
        outputs={},
        progress_timeout_hint=60,
    )


def _render(sc, target="s3_doublets_execute"):
    with tempfile.TemporaryDirectory() as tmp:
        run_dir = Path(tmp) / "run"
        repo_root = Path(tmp) / "repo"
        log = run_dir / "internal" / "hpc_monitor" / "logs" / "head.out"
        return render_submission_script(_spec(), sc, repo_root, run_dir, log, target)


class GpuRenderTests(unittest.TestCase):
    def test_cpu_default_exports_device_cpu_and_no_gpu_vars(self):
        sc = SiteConfig(scheduler="slurm", partition="cpu", account="vaquerizas",
                        conda_env="muagene")
        script = _render(sc)
        self.assertIn("export PMA_DEVICE=cpu", script)
        self.assertNotIn("PMA_SLURM_GPU_GRES", script)
        self.assertNotIn("PMA_CONDA_ENV_GPU", script)

    def test_slurm_gpu_vars_exported(self):
        sc = SiteConfig(scheduler="slurm", partition="cpu", account="vaquerizas",
                        conda_env="muagene", device="gpu", gpu_partition="gpu",
                        gpu_gres="gpu:A5000:1", gpu_conda_env="muagene-gpu")
        script = _render(sc)
        self.assertIn("export PMA_DEVICE=gpu", script)
        self.assertIn("export PMA_SLURM_GPU_PARTITION=gpu", script)
        self.assertIn("export PMA_SLURM_GPU_GRES=gpu:A5000:1", script)
        self.assertIn("export PMA_CONDA_ENV_GPU=muagene-gpu", script)

    def test_pbs_gpu_vars_exported(self):
        sc = SiteConfig(scheduler="pbs", queue="workq", project="vaquerizas",
                        conda_env="muagene", device="gpu",
                        gpu_select_extra="ngpus=1:gpu_type=a100", gpu_queue="gpuq",
                        gpu_conda_env="muagene-gpu")
        script = _render(sc)
        self.assertIn("export PMA_DEVICE=gpu", script)
        self.assertIn('export PMA_PBS_GPU_SELECT_EXTRA="ngpus=1:gpu_type=a100"', script)
        self.assertIn("export PMA_PBS_GPU_QUEUE=gpuq", script)
        self.assertIn("export PMA_CONDA_ENV_GPU=muagene-gpu", script)

    def test_load_site_config_round_trips_gpu_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "site.config"
            p.write_text(
                "schema_version: '1'\nscheduler: slurm\ndevice: gpu\n"
                "slurm:\n  partition: cpu\n  account: vaquerizas\n"
                "  gpu_partition: gpu\n  gpu_gres: gpu:A5000:1\n"
                "common:\n  conda_env: muagene\n  gpu_conda_env: muagene-gpu\n"
            )
            sc = load_site_config(p)
            self.assertEqual(sc.device, "gpu")
            self.assertEqual(sc.gpu_partition, "gpu")
            self.assertEqual(sc.gpu_gres, "gpu:A5000:1")
            self.assertEqual(sc.gpu_conda_env, "muagene-gpu")

    def test_legacy_site_config_without_device_defaults_cpu(self):
        # Backward-compat: existing site.configs have no device/gpu keys.
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "site.config"
            p.write_text(
                "schema_version: '1'\nscheduler: slurm\n"
                "slurm:\n  partition: cpu\n  account: vaquerizas\n"
                "common:\n  conda_env: grn\n"
            )
            sc = load_site_config(p)
            self.assertEqual(sc.device, "cpu")
            self.assertIsNone(sc.gpu_gres)
            self.assertIsNone(sc.gpu_conda_env)


if __name__ == "__main__":
    unittest.main()
