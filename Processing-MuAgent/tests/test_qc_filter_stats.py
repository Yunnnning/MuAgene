"""Tests for order-independent exclusive QC removal counts."""
from __future__ import annotations

import unittest

import numpy as np

from executor.methods.qc_filter_stats import append_frip_exclusive, exclusive_removals


class ExclusiveRemovalsTests(unittest.TestCase):
    def test_single_metric_failure(self):
        pass_a = np.array([True, True, False, True])
        pass_b = np.array([True, True, True, True])
        out = exclusive_removals({"a": pass_a, "b": pass_b})
        self.assertEqual(out["a"], 1)
        self.assertEqual(out["b"], 0)
        self.assertEqual(out["multiple_metrics"], 0)
        self.assertEqual(out["total_removed"], 1)

    def test_multiple_metric_failure(self):
        pass_a = np.array([True, False, False])
        pass_b = np.array([True, False, True])
        out = exclusive_removals({"a": pass_a, "b": pass_b})
        self.assertEqual(out["a"], 0)
        self.assertEqual(out["b"], 1)
        self.assertEqual(out["multiple_metrics"], 1)
        self.assertEqual(out["total_removed"], 2)

    def test_append_frip_updates_total(self):
        base = {"n_fragments": 2, "tss_enrichment": 1, "multiple_metrics": 0, "total_removed": 3}
        out = append_frip_exclusive(base, frip_fail=4, n_pre=100, n_post=93)
        self.assertEqual(out["frip_min"], 4)
        self.assertEqual(out["total_removed"], 7)
