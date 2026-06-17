"""Dry-run render tests for GPU routing in the SLURM/PBS profile submit scripts.

PBS has no hardware on this host, so its GPU directives are verified structurally
(render-equality via PMA_SUBMIT_DRY_RUN) exactly as the SLURM ones are. CPU paths
are asserted unchanged so existing runs are unaffected.
"""
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SLURM = REPO / "workflow" / "profiles" / "slurm" / "slurm-submit.sh"
PBS = REPO / "workflow" / "profiles" / "pbs" / "pbs-submit.sh"

# A run directory the GPU container wrapper must bind. Need not exist — DRY_RUN only
# prints the resolved wrapper text.
RUN_DIR = "/tmp/pma_test_run_dir"


def _run(script, gpu, env):
    """Run a submit script in dry-run mode; return the CompletedProcess."""
    with tempfile.TemporaryDirectory() as tmp:
        js = Path(tmp) / "snakejob.sh"
        js.write_text("#!/bin/sh\npython -m snakemake --mode subprocess x\n")
        full = {
            **os.environ,
            "PMA_SUBMIT_DRY_RUN": "1",
            "PMA_DISABLE_STORAGE_LOCAL_COPIES": "0",  # skip the python sanitize step
            "PMA_REPO_ROOT": str(REPO),
            "PMA_LOG_DIR": tmp,
            **env,
        }
        return subprocess.run(
            ["bash", str(script), "s3_doublets", "1", "2", "8000", "60", str(gpu), str(js)],
            capture_output=True, text=True, env=full,
        )


def _submit(script, gpu, env):
    r = _run(script, gpu, env)
    assert r.returncode == 0, r.stderr
    return r.stdout


class SlurmGpuRoutingTests(unittest.TestCase):
    def test_cpu_job_unchanged(self):
        out = _submit(SLURM, 0, {"PMA_SLURM_PARTITION": "cpu", "PMA_SLURM_ACCOUNT": "vaquerizas"})
        self.assertIn("--partition cpu", out)
        self.assertIn("--export=ALL", out)
        self.assertNotIn("--gres", out)
        self.assertNotIn("PMA_DEVICE=gpu", out)
        self.assertNotIn("--bind", out)  # no container wrapper for CPU jobs

    def test_gpu_conda_provider_routes_partition_gres_env(self):
        out = _submit(SLURM, 1, {
            "PMA_SLURM_PARTITION": "cpu", "PMA_SLURM_ACCOUNT": "vaquerizas",
            "PMA_SLURM_GPU_PARTITION": "gpu", "PMA_SLURM_GPU_GRES": "gpu:A5000:1",
            "PMA_CONDA_ENV_GPU": "muagene-gpu",
        })
        self.assertIn("--partition gpu", out)
        self.assertIn("--gres gpu:A5000:1", out)
        self.assertIn("PMA_DEVICE=gpu", out)
        self.assertIn("PMA_CONDA_ENV=muagene-gpu", out)

    def test_gpu_container_binds_repo_root_and_run_dir(self):
        out = _submit(SLURM, 1, {
            "PMA_SLURM_GPU_PARTITION": "gpu", "PMA_SLURM_GPU_GRES": "gpu:A5000:1",
            "PMA_GPU_PROVIDER": "container", "PMA_GPU_IMAGE": "/img/muagene-gpu.sif",
            "PMA_RUN_DIR": RUN_DIR,
        })
        self.assertIn("--gres gpu:A5000:1", out)
        self.assertIn("singularity exec --nv", out)
        self.assertIn("/img/muagene-gpu.sif", out)
        # The bug fix: the container must bind BOTH the repo root and the run dir.
        self.assertIn(f"--bind {REPO}", out)
        self.assertIn(f"--bind {RUN_DIR}", out)
        self.assertNotIn("PMA_CONDA_ENV=", out)  # container carries its own env

    def test_gpu_container_appends_optional_scratch_bind(self):
        out = _submit(SLURM, 1, {
            "PMA_SLURM_GPU_PARTITION": "gpu", "PMA_SLURM_GPU_GRES": "gpu:A5000:1",
            "PMA_GPU_PROVIDER": "container", "PMA_GPU_IMAGE": "/img/muagene-gpu.sif",
            "PMA_RUN_DIR": RUN_DIR, "PMA_GPU_BIND": "/scratch/fast",
        })
        self.assertIn(f"--bind {RUN_DIR}", out)
        self.assertIn("--bind /scratch/fast", out)

    def test_gpu_container_warns_when_run_dir_unset(self):
        r = _run(SLURM, 1, {
            "PMA_SLURM_GPU_PARTITION": "gpu", "PMA_SLURM_GPU_GRES": "gpu:A5000:1",
            "PMA_GPU_PROVIDER": "container", "PMA_GPU_IMAGE": "/img/muagene-gpu.sif",
            "PMA_RUN_DIR": "",  # force-empty so the script can't inherit a stray value
        })
        self.assertEqual(r.returncode, 0, r.stderr)         # warning, not a failure
        self.assertIn("PMA_RUN_DIR", r.stderr)
        self.assertIn("WARNING", r.stderr)
        self.assertIn(f"--bind {REPO}", r.stdout)            # repo still bound


class PbsGpuRoutingTests(unittest.TestCase):
    def test_cpu_job_unchanged(self):
        out = _submit(PBS, 0, {"PMA_PBS_QUEUE": "workq"})
        self.assertIn("select=1:ncpus=2:mem=8000mb", out)
        self.assertNotIn("ngpus", out)
        self.assertNotIn("PMA_DEVICE=gpu", out)
        self.assertNotIn("--bind", out)

    def test_gpu_appends_ngpus_and_routes_queue(self):
        out = _submit(PBS, 1, {
            "PMA_PBS_QUEUE": "workq", "PMA_PBS_GPU_SELECT_EXTRA": "ngpus=1:gpu_type=a100",
            "PMA_PBS_GPU_QUEUE": "gpuq", "PMA_CONDA_ENV_GPU": "muagene-gpu",
        })
        self.assertIn("select=1:ncpus=2:mem=8000mb:ngpus=1:gpu_type=a100", out)
        self.assertIn("-q gpuq", out)
        self.assertIn("PMA_DEVICE=gpu", out)
        self.assertIn("PMA_CONDA_ENV=muagene-gpu", out)

    def test_gpu_container_binds_repo_root_and_run_dir(self):
        out = _submit(PBS, 1, {
            "PMA_PBS_GPU_SELECT_EXTRA": "ngpus=1", "PMA_GPU_PROVIDER": "container",
            "PMA_GPU_IMAGE": "/img/muagene-gpu.sif", "PMA_RUN_DIR": RUN_DIR,
        })
        self.assertIn(":ngpus=1", out)
        self.assertIn("singularity exec --nv", out)
        self.assertIn("/img/muagene-gpu.sif", out)
        self.assertIn(f"--bind {REPO}", out)
        self.assertIn(f"--bind {RUN_DIR}", out)

    def test_gpu_container_warns_when_run_dir_unset(self):
        r = _run(PBS, 1, {
            "PMA_PBS_GPU_SELECT_EXTRA": "ngpus=1", "PMA_GPU_PROVIDER": "container",
            "PMA_GPU_IMAGE": "/img/muagene-gpu.sif", "PMA_RUN_DIR": "",
        })
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("PMA_RUN_DIR", r.stderr)
        self.assertIn("WARNING", r.stderr)


if __name__ == "__main__":
    unittest.main()
