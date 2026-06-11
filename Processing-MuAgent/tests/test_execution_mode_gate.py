"""Tests for the mandatory execution-mode confirmation gate.

System requirement: Processing-MuAgent must always confirm execution mode
(local vs HPC) with the user before launching ANY compute job — enforced
unconditionally on every `run` and `submit`, including resume submissions and
runs whose config never recorded an execution mode.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import yaml
from click.testing import CliRunner

from executor import cli
from executor import provenance as _prov
from executor.run_paths import RunPaths


def _init_run(tmp: Path, *, mode: str | None = None, confirmed: bool = False) -> Path:
    """Scaffold a run dir via `init`; optionally seed execution mode/confirmation.

    Returns the canonical config path (deliverables/plan/config/run.yaml).
    """
    run_dir = tmp / "run"
    draft = tmp / "run.yaml.draft"
    draft.write_text(yaml.safe_dump({"run_dir": str(run_dir)}))
    res = CliRunner().invoke(cli.main, ["init", "--config", str(draft)])
    assert res.exit_code == 0, res.output
    paths = RunPaths(run_dir)
    if mode is not None:
        _prov.set_param(str(paths.parameters_yaml), "execution.mode", mode,
                        source="user", confidence="high", rationale="test")
    if confirmed:
        _prov.set_param(str(paths.parameters_yaml), "execution.user_confirmed", True,
                        source="user", confidence="high", rationale="test")
    return paths.run_yaml


class GateHelperTests(unittest.TestCase):
    def test_raises_when_mode_unset(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _init_run(Path(tmp))
            run_dir = cfg.parents[3]
            with self.assertRaises(cli.click.ClickException) as ctx:
                cli._enforce_execution_mode_gate(run_dir, RunPaths(run_dir))
            self.assertIn("not set", ctx.exception.message)
            self.assertIn("configure-execution", ctx.exception.message)

    def test_raises_when_set_but_not_confirmed(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _init_run(Path(tmp), mode="local")
            run_dir = cfg.parents[3]
            with self.assertRaises(cli.click.ClickException) as ctx:
                cli._enforce_execution_mode_gate(run_dir, RunPaths(run_dir))
            self.assertIn("not confirmed", ctx.exception.message)

    def test_passes_when_confirmed(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _init_run(Path(tmp), mode="local", confirmed=True)
            run_dir = cfg.parents[3]
            cli._enforce_execution_mode_gate(run_dir, RunPaths(run_dir))  # no raise

    def test_idempotent_on_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _init_run(Path(tmp), mode="slurm", confirmed=True)
            run_dir = cfg.parents[3]
            cli._enforce_execution_mode_gate(run_dir, RunPaths(run_dir))
            cli._enforce_execution_mode_gate(run_dir, RunPaths(run_dir))  # still no raise


class ConfigureExecutionFlagTests(unittest.TestCase):
    def test_confirmed_by_user_records_true(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _init_run(Path(tmp))
            res = CliRunner().invoke(cli.main, [
                "configure-execution", "--config", str(cfg),
                "--mode", "local", "--confirmed-by-user"])
            self.assertEqual(res.exit_code, 0, res.output)
            params = RunPaths(cfg.parents[3]).parameters_yaml
            self.assertIs(_prov.get_value(str(params), "execution.user_confirmed"), True)

    def test_default_not_confirmed_false_and_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _init_run(Path(tmp))
            res = CliRunner().invoke(cli.main, [
                "configure-execution", "--config", str(cfg), "--mode", "local"])
            self.assertEqual(res.exit_code, 0, res.output)
            params = RunPaths(cfg.parents[3]).parameters_yaml
            self.assertIs(_prov.get_value(str(params), "execution.user_confirmed"), False)
            self.assertIn("NOT user-confirmed", res.output)

    def test_reconfig_same_mode_preserves_confirmation(self):
        # Confirm once, then re-config the SAME mode without the flag (e.g. a
        # resource tweak) — confirmation must be preserved, not dropped.
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _init_run(Path(tmp))
            run = CliRunner()
            run.invoke(cli.main, ["configure-execution", "--config", str(cfg),
                                  "--mode", "local", "--confirmed-by-user"])
            res = run.invoke(cli.main, ["configure-execution", "--config", str(cfg),
                                        "--mode", "local"])
            self.assertEqual(res.exit_code, 0, res.output)
            params = RunPaths(cfg.parents[3]).parameters_yaml
            self.assertIs(_prov.get_value(str(params), "execution.user_confirmed"), True)

    def test_reconfig_changed_mode_resets_confirmation(self):
        # Confirm local, then switch to slurm without the flag — a mode change
        # must reset to unconfirmed so the user re-confirms the new backend.
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _init_run(Path(tmp))
            run = CliRunner()
            run.invoke(cli.main, ["configure-execution", "--config", str(cfg),
                                  "--mode", "local", "--confirmed-by-user"])
            res = run.invoke(cli.main, ["configure-execution", "--config", str(cfg),
                                        "--mode", "slurm", "--slurm-partition", "cpu"])
            self.assertEqual(res.exit_code, 0, res.output)
            params = RunPaths(cfg.parents[3]).parameters_yaml
            self.assertIs(_prov.get_value(str(params), "execution.user_confirmed"), False)


class RunGateCliTests(unittest.TestCase):
    @mock.patch("executor.cli._snakemake")
    def test_run_blocks_before_snakemake_when_unconfirmed(self, m_snake):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _init_run(Path(tmp), mode="local")  # not confirmed
            res = CliRunner().invoke(cli.main, [
                "run", "--config", str(cfg), "--target", "s0_ingest_execute",
                "--no-context"])
            self.assertNotEqual(res.exit_code, 0)
            self.assertIn("not confirmed", res.output)
            m_snake.assert_not_called()

    @mock.patch("executor.cli._snakemake")
    def test_run_refuses_cluster_mode_local_only(self, m_snake):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _init_run(Path(tmp), mode="pbs", confirmed=True)
            res = CliRunner().invoke(cli.main, [
                "run", "--config", str(cfg), "--target", "s0_ingest_execute",
                "--no-context"])
            self.assertNotEqual(res.exit_code, 0)
            self.assertIn("local-only", res.output)
            self.assertIn("submit", res.output)
            self.assertIn("--executor pbs", res.output)
            m_snake.assert_not_called()


class SubmitGateCliTests(unittest.TestCase):
    def test_submit_blocks_when_mode_unset(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _init_run(Path(tmp))  # no mode at all
            res = CliRunner().invoke(cli.main, [
                "submit", "--config", str(cfg), "--executor", "slurm",
                "--target", "s0_ingest_execute"])
            self.assertNotEqual(res.exit_code, 0)
            self.assertIn("not set", res.output)


if __name__ == "__main__":
    unittest.main()
