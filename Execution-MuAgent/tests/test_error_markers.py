"""Error markers match whole tokens: a false workflow_error_marker kills a healthy job."""
from __future__ import annotations

import unittest

from execution_muagent.monitor import _extract_markers


class ErrorMarkerMatchingTests(unittest.TestCase):
    def test_words_containing_oom_do_not_match(self):
        benign = "Saved figure to room_plan.png; zoom level 3; bloom filter built; Doom.\n"
        self.assertEqual(_extract_markers(benign), [])

    def test_real_oom_signals_still_match(self):
        for line in (
            "slurmstepd: error: Detected 1 oom_kill event in StepId=123.batch",
            "slurmstepd: error: Detected 1 oom-kill event(s)",
            "Job 123 OOM: memory limit exceeded",
            "State: OUT_OF_MEMORY",
            "RuntimeError: CUDA out of memory.",
        ):
            self.assertTrue(_extract_markers(line), line)

    def test_other_markers_unchanged(self):
        log = ("Error in rule s1_rna_qc_execute:\nTraceback (most recent call last):\n"
               "FileNotFoundError: x\nslurmstepd: *** JOB 9 CANCELLED AT 2026 DUE TO TIME LIMIT ***")
        found = set(_extract_markers(log))
        self.assertTrue({"Error in rule", "Traceback", "FileNotFoundError",
                         "CANCELLED", "DUE TO TIME LIMIT"} <= found)


if __name__ == "__main__":
    unittest.main()
