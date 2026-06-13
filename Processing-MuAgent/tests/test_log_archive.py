"""Archiving the prior run's Snakemake logs on resubmit (executor/hpc.py).

Covers the stale-state bug where, right after a resubmit, `hpc-status` reported a
phantom failure because `collect_failed_snakemake_rules` read the *previous* run's
per-rule + main logs while the new head job was still PENDING. `archive_prior_run_logs`
moves those logs aside so the live dirs only ever describe the current run; the fix
preserves history (move, not delete) and leaves Snakemake's own state intact.
"""
import tempfile
import unittest
from pathlib import Path

from executor import hpc
from executor.run_paths import RunPaths
from executor import stage_progress


def _seed_failed_logs(sm_root: Path, rule: str = "s5_atac_spectral_execute") -> None:
    (sm_root / "log").mkdir(parents=True, exist_ok=True)
    (sm_root / "log" / "2026-06-13T1.snakemake.log").write_text(
        "Error in rule s5_atac_spectral_execute:\nExiting because a job execution failed\n"
    )
    rd = sm_root / "slurm_logs" / f"rule_{rule}"
    rd.mkdir(parents=True, exist_ok=True)
    (rd / "999.log").write_text("RuleException\nPanicException\n")
    (sm_root / "metadata").mkdir(exist_ok=True)  # snakemake state — must survive


class ArchivePriorRunLogsTests(unittest.TestCase):
    def test_moves_logs_and_preserves_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            wd = Path(tmp) / "snakemake"
            sm = wd / ".snakemake"
            _seed_failed_logs(sm)
            arch = hpc.archive_prior_run_logs(wd, stamp="20260613T160000Z")
            self.assertIsNotNone(arch)
            self.assertTrue((arch / "log" / "2026-06-13T1.snakemake.log").exists())
            self.assertTrue((arch / "slurm_logs" / "rule_s5_atac_spectral_execute" / "999.log").exists())
            self.assertEqual(list((sm / "log").glob("*.log")), [])
            self.assertFalse((sm / "slurm_logs").exists())
            self.assertTrue((sm / "metadata").exists())  # untouched

    def test_noop_when_nothing_to_archive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            wd = Path(tmp) / "snakemake"
            (wd / ".snakemake").mkdir(parents=True)
            self.assertIsNone(hpc.archive_prior_run_logs(wd))

    def test_archived_failure_no_longer_marks_stage_failed(self) -> None:
        # End-to-end: a prior failure is invisible to stage-state derivation
        # once archived (the new run's PENDING window shows no phantom failure).
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            paths = RunPaths(run_dir)
            sm = paths.snakemake_workdir / ".snakemake"
            _seed_failed_logs(sm)
            self.assertIn(
                "s5_atac_spectral_execute",
                stage_progress.collect_failed_snakemake_rules(paths),
            )
            hpc.archive_prior_run_logs(paths.snakemake_workdir, stamp="20260613T160000Z")
            self.assertEqual(
                stage_progress.collect_failed_snakemake_rules(paths), frozenset()
            )


if __name__ == "__main__":
    unittest.main()
