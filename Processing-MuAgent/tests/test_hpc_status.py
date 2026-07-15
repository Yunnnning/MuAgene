"""Tests for report-and-repoll surfacing in `hpc-status` and head-job spec outputs."""
from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path

import yaml
from click.testing import CliRunner

from executor import cli, specs
from executor.run_paths import RunPaths


def _setup_run(tmp: str) -> tuple[Path, Path]:
    run_dir = Path(tmp) / "run"
    paths = RunPaths(run_dir)
    paths.ensure()
    cfg = run_dir / "run.yaml"
    cfg.write_text(yaml.safe_dump({"run_dir": str(run_dir)}), encoding="utf-8")
    mon = run_dir / "internal" / "hpc_monitor"
    mon.mkdir(parents=True, exist_ok=True)
    (mon / "latest_submission.json").write_text(json.dumps({
        "job_id": "1015707", "target": "s0_ingest_execute",
        "submitted_at": "2026-06-12T13:12:44Z",
    }), encoding="utf-8")
    return run_dir, cfg


def _write_snapshot(run_dir: Path, state: str) -> None:
    mon = run_dir / "internal" / "hpc_monitor"
    (mon / "latest_snapshot.json").write_text(json.dumps({
        "scheduler": {"state": state, "elapsed": "00:04:22", "timelimit": "1-00:00:00"},
        "monitor_state": {"health": "healthy", "silence_intervals": 0, "tolerance_n": 27},
        "findings": [],
        "interval_s": 270.0,
        "next_recheck_after_s": 295.0,
    }), encoding="utf-8")


class HpcStatusReportAndRepollTests(unittest.TestCase):
    def test_running_job_prints_next_check_and_fingerprint(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir, cfg = _setup_run(tmp)
            _write_snapshot(run_dir, "RUNNING")
            # A live monitor.pid (this process) → supervisor "alive" → re-poll path.
            (run_dir / "internal" / "hpc_monitor" / "monitor.pid").write_text(str(os.getpid()))
            res = CliRunner().invoke(cli.main, ["hpc-status", "--config", str(cfg)])
            self.assertEqual(res.exit_code, 0, res.output)
            self.assertIn("Next check: re-poll", res.output)
            self.assertIn("~295s", res.output)
            self.assertIn("State: RUNNING/healthy/sup=alive/gate=none/findings=0", res.output)
            self.assertNotIn("Gate signal present", res.output)

    def test_gate_armed_prints_gate_signal_and_no_repoll(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir, cfg = _setup_run(tmp)
            _write_snapshot(run_dir, "COMPLETED")
            # No monitor.pid (daemon exited) + a review gate armed.
            RunPaths(run_dir).awaiting_sentinel("plan_review").write_text("")
            res = CliRunner().invoke(cli.main, ["hpc-status", "--config", str(cfg)])
            self.assertEqual(res.exit_code, 0, res.output)
            self.assertIn("Gate signal present", res.output)
            self.assertIn("gate=awaiting_approval", res.output)
            self.assertNotIn("Next check: re-poll", res.output)


class HeadJobSpecOutputsTests(unittest.TestCase):
    def test_s0_execute_target_populates_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = specs.write_head_job_spec(tmp, "s0_ingest_execute")
            spec = yaml.safe_load(p.read_text())
            # Head-job outputs are the durable declared markers the monitor can verify
            # (validation_report.json + plan), never the deletable rna_ingest.h5ad cache
            # (which is not a declared Snakemake output and is removed at the QC gate).
            self.assertIn("validation_report", spec["outputs"])
            self.assertIn("plan", spec["outputs"])
            self.assertNotIn("rna_h5ad", spec["outputs"])

    def test_multistage_target_leaves_outputs_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            for tgt in ("post_qc_review_propose", "all"):
                p = specs.write_head_job_spec(tmp, tgt)
                spec = yaml.safe_load(p.read_text())
                self.assertEqual(spec["outputs"], {}, tgt)


if __name__ == "__main__":
    unittest.main()
