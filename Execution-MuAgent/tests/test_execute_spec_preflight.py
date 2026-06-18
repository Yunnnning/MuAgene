"""execute-spec env preflight must fail loud.

A reconcile that RAISES (corrupt env_state.json, probe/subprocess crash) or RETURNS
error findings must abort the submit — never silently degrade to a job running against
an unverified env. submit_from_spec must not be reached in either case.
"""
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from click.testing import CliRunner

from execution_muagent import cli
from execution_muagent.monitor import SiteConfig


def _invoke(reconcile_side_effect=None, reconcile_return=None):
    """Run `execute-spec` with the spec/site loaders + submit stubbed out.

    Returns (result, sentinel) where sentinel['submitted'] records whether
    submit_from_spec was reached.
    """
    sentinel = {"submitted": False}

    def _fake_submit(*args, **kwargs):
        sentinel["submitted"] = True
        return {"rejected_as": None, "stdout": "", "stderr": "", "job_id": "1"}

    rec = mock.Mock()
    if reconcile_side_effect is not None:
        rec.side_effect = reconcile_side_effect
    else:
        rec.return_value = reconcile_return

    with tempfile.TemporaryDirectory() as tmp:
        spec = Path(tmp) / "head_job.yaml"
        spec.write_text("stage: s_test\n")
        site = Path(tmp) / "site.config"
        site.write_text("scheduler: slurm\n")
        repo = Path(tmp) / "repo"
        repo.mkdir()
        run_dir = Path(tmp) / "run"

        site_cfg = SiteConfig(scheduler="slurm", partition="cpu",
                              account="vaquerizas", conda_env="muagene")
        runner = CliRunner()
        with mock.patch.object(cli, "load_stage_spec",
                               return_value=mock.Mock(stage="s_test")), \
             mock.patch.object(cli, "load_site_config", return_value=site_cfg), \
             mock.patch.object(cli, "validate_spec", return_value=[]), \
             mock.patch("execution_muagent.environment.reconcile", rec), \
             mock.patch.object(cli, "submit_from_spec", side_effect=_fake_submit):
            result = runner.invoke(cli.main, [
                "execute-spec", "--spec", str(spec), "--site-config", str(site),
                "--run-dir", str(run_dir), "--repo-root", str(repo),
            ])
    return result, sentinel


class ExecuteSpecPreflightTests(unittest.TestCase):
    def test_crashed_reconcile_aborts_submit(self):
        result, sentinel = _invoke(reconcile_side_effect=RuntimeError("probe blew up"))
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("preflight", result.output.lower())
        self.assertFalse(sentinel["submitted"],
                         "submit must not run after a crashed reconcile")

    def test_reconcile_error_findings_abort_submit(self):
        result, sentinel = _invoke(reconcile_return={
            "ok": False,
            "findings": [{"severity": "error", "message": "cpu env is missing"}],
        })
        self.assertNotEqual(result.exit_code, 0)
        self.assertFalse(sentinel["submitted"],
                         "submit must not run when reconcile returns error findings")


def _devices_reconciled(device, gpu_stages_present):
    """Return the device list execute-spec actually reconciles, given a site.config
    device and a spec.gpu_stages_present. reconcile is stubbed to succeed; submit is
    stubbed to stop right after preflight so only the device gating is exercised."""
    rec = mock.Mock(return_value={"ok": True, "findings": [], "provision": None})
    with tempfile.TemporaryDirectory() as tmp:
        spec_path = Path(tmp) / "head_job.yaml"
        spec_path.write_text("stage: s_test\n")
        site = Path(tmp) / "site.config"
        site.write_text("scheduler: slurm\n")
        repo = Path(tmp) / "repo"
        repo.mkdir()
        run_dir = Path(tmp) / "run"
        site_cfg = SiteConfig(scheduler="slurm", partition="cpu", account="vaquerizas",
                              conda_env="muagene", device=device)
        spec = mock.Mock(stage="s_test", gpu_stages_present=gpu_stages_present)
        runner = CliRunner()
        with mock.patch.object(cli, "load_stage_spec", return_value=spec), \
             mock.patch.object(cli, "load_site_config", return_value=site_cfg), \
             mock.patch.object(cli, "validate_spec", return_value=[]), \
             mock.patch("execution_muagent.environment.reconcile", rec), \
             mock.patch.object(cli, "submit_from_spec",
                               side_effect=RuntimeError("stop after preflight")):
            runner.invoke(cli.main, [
                "execute-spec", "--spec", str(spec_path), "--site-config", str(site),
                "--run-dir", str(run_dir), "--repo-root", str(repo),
            ])
    return [c.args[2] for c in rec.call_args_list]


class ExecuteSpecGpuGatingTests(unittest.TestCase):
    def test_gpu_device_without_gpu_stages_reconciles_cpu_only(self):
        # device=gpu but no GPU-capable stage in the run → no GPU env pull.
        self.assertEqual(_devices_reconciled("gpu", False), ["cpu"])

    def test_gpu_device_with_gpu_stages_reconciles_both(self):
        self.assertEqual(_devices_reconciled("gpu", True), ["cpu", "gpu"])

    def test_cpu_device_reconciles_cpu_only(self):
        self.assertEqual(_devices_reconciled("cpu", True), ["cpu"])


if __name__ == "__main__":
    unittest.main()
