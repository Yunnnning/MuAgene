"""Tests for _seed_approvals idempotency.

A `submit --auto-approve` re-enters phases whose gates are already approved. It
must NOT re-stamp an existing <gate>.approved sentinel: rewriting it bumps the
file mtime, and because every QC/downstream execute rule declares the sentinel as
an input, a fresh mtime is an "Updated input files" trigger that needlessly
re-runs the already-approved upstream chain.
"""
from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

import yaml

from executor import approval, cli
from executor.run_paths import RunPaths


def _init_run(tmp: str) -> RunPaths:
    paths = RunPaths(tmp)
    paths.ensure()
    paths.parameters_yaml.write_text(
        yaml.safe_dump({"plan": {"workflow_branch": "paired"}})
    )
    return paths


class SeedApprovalsIdempotencyTests(unittest.TestCase):
    def setUp(self):
        self._had_flag = "PMA_AUTO_APPROVE" in os.environ
        self._prev_flag = os.environ.get("PMA_AUTO_APPROVE")

    def tearDown(self):
        if self._had_flag:
            os.environ["PMA_AUTO_APPROVE"] = self._prev_flag  # type: ignore[assignment]
        else:
            os.environ.pop("PMA_AUTO_APPROVE", None)

    def test_already_approved_sentinel_not_restamped(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            _init_run(tmp)
            approval.approve(run_dir, "post_qc_review", note="real approval")
            sentinel = RunPaths(run_dir).approved_sentinel("post_qc_review")
            before_bytes = sentinel.read_bytes()
            before_mtime = sentinel.stat().st_mtime_ns

            seeded = cli._seed_approvals(
                run_dir, ("post_qc_review",), note="auto-approved (submit)")

            self.assertEqual(seeded, [], "an already-approved gate must not be re-seeded")
            self.assertEqual(sentinel.read_bytes(), before_bytes,
                             "sentinel content (incl. original note/timestamp) must be preserved")
            self.assertEqual(sentinel.stat().st_mtime_ns, before_mtime,
                             "sentinel mtime must not be bumped — that would re-run the QC chain")

    def test_unapproved_gate_is_seeded(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            _init_run(tmp)
            seeded = cli._seed_approvals(
                run_dir, ("s7_clustering",), note="auto-approved (submit)")
            self.assertEqual(seeded, ["s7_clustering"])
            self.assertTrue(approval.is_approved(run_dir, "s7_clustering"))

    def test_protection_flag_set_even_when_nothing_freshly_seeded(self):
        # All gates already approved → nothing seeded, but PMA_AUTO_APPROVE must
        # still be set so the head-job's propose rules don't revoke the approvals.
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            _init_run(tmp)
            approval.approve(run_dir, "post_qc_review", note="real")
            os.environ.pop("PMA_AUTO_APPROVE", None)

            seeded = cli._seed_approvals(
                run_dir, ("post_qc_review",), note="auto-approved (submit)")

            self.assertEqual(seeded, [])
            self.assertEqual(os.environ.get("PMA_AUTO_APPROVE"), "1")

    def test_kept_gate_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            _init_run(tmp)
            seeded = cli._seed_approvals(
                run_dir, ("post_qc_review", "s7_clustering"),
                note="auto-approved (submit)", kept={"s7_clustering"})
            self.assertEqual(seeded, ["post_qc_review"])
            self.assertFalse(approval.is_approved(run_dir, "s7_clustering"))


if __name__ == "__main__":
    unittest.main()
