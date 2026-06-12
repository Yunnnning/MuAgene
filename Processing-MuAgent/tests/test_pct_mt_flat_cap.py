"""Regression test: the RNA pct_mt ceiling is a single flat cap (20%).

The sample-type-specific cap (20% cells / 10% nuclei) was removed because the
cap is only a guardrail bounding the MAD-derived threshold, and cells-vs-nuclei
is often ambiguous in user-declared context. The cap must now be 20% for every
sample type.
"""
from __future__ import annotations

import unittest

from executor import plan_assembler as pa


class PctMtFlatCapTests(unittest.TestCase):
    def test_flat_20_cap_for_all_sample_types(self):
        for sample_type in ("cells", "nuclei", "unknown", "tissue"):
            plan = pa.assemble_plan(
                "/tmp/_pct_mt_flat_cap_test",
                workflow_branch="rna_only",
                sample_type=sample_type,
            )
            ceil = plan["stages"]["s1_rna_qc"]["parameters"]["pct_mt_ceiling"]
            self.assertEqual(ceil["value"], 20.0,
                             f"pct_mt cap must be 20 for sample_type={sample_type}")
            self.assertNotIn("nuclei", ceil["rationale"].lower(),
                             "cap rationale must not be sample-type-specific")


if __name__ == "__main__":
    unittest.main()
