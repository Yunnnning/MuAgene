"""Tests for plan-review appendix stage ordering."""
from __future__ import annotations

import unittest

from executor.plan_assembler import render_plan_appendix


class PlanAppendixOrderTests(unittest.TestCase):
    def test_s1a_before_s1_rna_qc_in_appendix(self):
        plan = {
            "workflow_branch": "paired",
            "stages": {
                "s1_rna_qc": {"parameters": {"total_counts_k_mad": {"value": 5.0, "rationale": "x"}}},
                "s1a_ambient": {"parameters": {"method": {"value": "auto", "rationale": "y"}}},
            },
            "warnings": [],
        }
        md = render_plan_appendix(plan)
        s1a_pos = md.index("### s1a_ambient")
        s1_pos = md.index("### s1_rna_qc")
        self.assertLess(s1a_pos, s1_pos)


if __name__ == "__main__":
    unittest.main()
