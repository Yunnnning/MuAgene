"""Tests for order-independent marginal QC removal counts."""
from __future__ import annotations

import unittest

import numpy as np

from executor.methods.qc_filter_stats import append_frip, marginal_removals


class MarginalRemovalsTests(unittest.TestCase):
    def test_single_metric_failure(self):
        pass_a = np.array([True, True, False, True])
        pass_b = np.array([True, True, True, True])
        out = marginal_removals({"a": pass_a, "b": pass_b})
        self.assertEqual(out["a"], 1)
        self.assertEqual(out["b"], 0)
        self.assertEqual(out["multiple_metrics"], 0)
        self.assertEqual(out["total_removed"], 1)

    def test_marginal_counts_each_failure(self):
        # idx0 passes both; idx1 fails both; idx2 fails only a.
        pass_a = np.array([True, False, False])
        pass_b = np.array([True, False, True])
        out = marginal_removals({"a": pass_a, "b": pass_b})
        # Marginal: every cell failing the threshold is counted, so idx1 + idx2.
        self.assertEqual(out["a"], 2)
        self.assertEqual(out["b"], 1)   # idx1 fails b
        self.assertEqual(out["multiple_metrics"], 1)  # idx1 fails >= 2
        self.assertEqual(out["total_removed"], 2)     # union: idx1, idx2

    def test_order_independent(self):
        masks = {
            "a": np.array([True, False, False, True]),
            "b": np.array([False, False, True, True]),
        }
        out1 = marginal_removals(masks)
        out2 = marginal_removals({"b": masks["b"], "a": masks["a"]})
        self.assertEqual(out1["a"], out2["a"])
        self.assertEqual(out1["b"], out2["b"])
        self.assertEqual(out1["multiple_metrics"], out2["multiple_metrics"])
        self.assertEqual(out1["total_removed"], out2["total_removed"])

    def test_empty(self):
        out = marginal_removals({})
        self.assertEqual(out["multiple_metrics"], 0)
        self.assertEqual(out["total_removed"], 0)

    def test_append_frip_updates_total(self):
        base = {"n_fragments": 2, "tss_enrichment": 1, "multiple_metrics": 0, "total_removed": 3}
        out = append_frip(base, frip_fail=4, n_pre=100, n_post=93)
        self.assertEqual(out["frip_min"], 4)
        self.assertEqual(out["total_removed"], 7)

    def test_append_frip_none_keeps_union(self):
        base = {"n_fragments": 2, "multiple_metrics": 0, "total_removed": 3}
        out = append_frip(base, frip_fail=None, n_pre=100, n_post=97)
        self.assertNotIn("frip_min", out)
        self.assertEqual(out["total_removed"], 3)


if __name__ == "__main__":
    unittest.main()
