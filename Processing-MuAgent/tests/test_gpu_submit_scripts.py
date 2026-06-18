"""Dry-run render tests for the SLURM/PBS profile submit scripts.

Preprocessing is CPU-only (_GPU_CAPABLE is empty in resources.smk). The submit
scripts no longer accept a gpu argument or perform GPU routing — that belongs in the
integration pipeline's submit profile. These tests verify:

  - CPU jobs are submitted to the expected partition with correct resources.
  - PMA_DEVICE=cpu is always exported, overriding any GPU value carried by the
    head-job (head-job may be configured --device gpu for future integration).
"""
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SLURM = REPO / "workflow" / "profiles" / "slurm" / "slurm-submit.sh"
PBS = REPO / "workflow" / "profiles" / "pbs" / "pbs-submit.sh"


def _run(script, env):
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
            ["bash", str(script), "s0_ingest_execute", "1", "4", "128000", "360", str(js)],
            capture_output=True, text=True, env=full,
        )


def _submit(script, env):
    r = _run(script, env)
    assert r.returncode == 0, r.stderr
    return r.stdout


class SlurmSubmitTests(unittest.TestCase):
    def test_cpu_job_routed_correctly(self):
        out = _submit(SLURM, {"PMA_SLURM_PARTITION": "cpu", "PMA_SLURM_ACCOUNT": "vaquerizas"})
        self.assertIn("--partition cpu", out)
        self.assertIn("--account vaquerizas", out)
        self.assertIn("--cpus-per-task=4", out)
        self.assertIn("--mem=128000M", out)
        self.assertNotIn("--gres", out)

    def test_pma_device_cpu_always_exported(self):
        # Even when head-job has PMA_DEVICE=gpu, child preprocessing jobs must
        # always get PMA_DEVICE=cpu (preprocessing is CPU-only).
        out = _submit(SLURM, {
            "PMA_SLURM_PARTITION": "cpu", "PMA_SLURM_ACCOUNT": "vaquerizas",
            "PMA_DEVICE": "gpu",
        })
        self.assertIn("PMA_DEVICE=cpu", out)
        self.assertNotIn("PMA_DEVICE=gpu", out)


class PbsSubmitTests(unittest.TestCase):
    def test_cpu_job_routed_correctly(self):
        out = _submit(PBS, {"PMA_PBS_QUEUE": "workq", "PMA_PBS_PROJECT": "vaquerizas"})
        self.assertIn("select=1:ncpus=4:mem=128000mb", out)
        self.assertNotIn("ngpus", out)
        self.assertIn("-q workq", out)
        self.assertIn("-P vaquerizas", out)

    def test_pma_device_cpu_always_exported(self):
        out = _submit(PBS, {"PMA_PBS_QUEUE": "workq", "PMA_DEVICE": "gpu"})
        self.assertIn("PMA_DEVICE=cpu", out)
        self.assertNotIn("PMA_DEVICE=gpu", out)


if __name__ == "__main__":
    unittest.main()
