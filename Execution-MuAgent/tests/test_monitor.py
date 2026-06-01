"""Unit tests for HPC monitor helpers (no scratch dirs in the package tree)."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from execution_muagent.monitor import (
    Submission,
    _is_run_scoped_progress_path,
    discover_child_job_ids,
    parse_job_ids_from_log,
)


class MonitorHelperTests(unittest.TestCase):
    def test_parse_job_ids_snakemake_external_jobid(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "snakemake.log"
            log.write_text(
                "Submitted job 10 with external jobid '992575'.\n"
                "Submitted job 2 with external jobid '992578'.\n",
                encoding="utf-8",
            )
            self.assertEqual(parse_job_ids_from_log(log), ["992575", "992578"])

    def test_parse_job_ids_from_slurm_log_filename(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "992579.log"
            log.write_text("Finished jobid: 0\n", encoding="utf-8")
            self.assertEqual(parse_job_ids_from_log(log), ["992579"])

    def test_discover_child_job_ids_from_run_logs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_dir = root / "run"
            slurm_logs = (
                run_dir / "internal" / "snakemake" / ".snakemake" / "slurm_logs" / "rule_x"
            )
            slurm_logs.mkdir(parents=True)
            (slurm_logs / "992578.log").write_text(
                "Submitted job 2 with external jobid '992578'.\n",
                encoding="utf-8",
            )
            submission = Submission(
                agent="Processing-MuAgent",
                executor="slurm",
                job_id="992574",
                run_dir=str(run_dir),
                config=str(run_dir / "cfg.yaml"),
                target="all",
                repo_root=str(root / "repo"),
                log_path=str(root / "repo" / "logs" / "pma_runner-992574.out"),
                submitted_at="2026-06-01T11:57:46Z",
            )
            self.assertEqual(discover_child_job_ids(submission), ["992578"])

    def test_run_scoped_progress_excludes_other_run_head_log(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = root / "run_a"
            run_b = root / "run_b"
            run_a.mkdir()
            run_b.mkdir()
            other_log = root / "logs" / "pma_runner-999.out"
            other_log.parent.mkdir()
            other_log.write_text("other run\n", encoding="utf-8")
            submission = Submission(
                agent="Processing-MuAgent",
                executor="slurm",
                job_id="992574",
                run_dir=str(run_a),
                config=str(run_a / "cfg.yaml"),
                target="all",
                repo_root=str(root / "repo"),
                log_path=str(root / "logs" / "pma_runner-992574.out"),
                submitted_at="2026-06-01T11:57:46Z",
            )
            self.assertFalse(_is_run_scoped_progress_path(submission, other_log))
            self.assertTrue(
                _is_run_scoped_progress_path(submission, run_a / "internal" / "log.jsonl")
            )


if __name__ == "__main__":
    unittest.main()
