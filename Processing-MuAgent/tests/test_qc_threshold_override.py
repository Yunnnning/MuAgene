"""Tests for the manual QC threshold-override layer.

Overrides pin the EFFECTIVE filtering bound to any value the user chooses, while
the MAD/floor derivation still runs and is exposed as the ``*_derived`` keys (grey
reference lines). The chosen override is drawn red. The no-override path must stay
byte-identical to the previous behaviour.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd

from executor import figures
from executor import provenance as _prov
from executor.methods import qc_thresholds as qct


def _rna_obs(n: int = 500, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "total_counts": rng.integers(50, 5000, n).astype(float),
        "n_genes_by_counts": rng.integers(20, 2000, n).astype(float),
        "pct_counts_mt": rng.uniform(0, 8, n),
        "pct_counts_ribo": rng.uniform(0, 30, n),
    })


_RNA_KW = dict(
    total_counts_k_mad=5.0, n_genes_k_mad=5.0, pct_mt_k=3.0,
    pct_mt_ceiling=20.0, pct_mt_floor=5.0,
    min_counts_floor=500, min_genes_floor=250,
)


class RnaOverrideMathTests(unittest.TestCase):
    def test_override_sets_effective_bound_and_keeps_derived(self):
        obs = _rna_obs()
        base = qct.rna_thresholds(obs, **_RNA_KW)
        th = qct.rna_thresholds(obs, **_RNA_KW, n_genes_min_override=300)
        # Effective bound is exactly the override; derived is the original MAD/floor.
        self.assertEqual(th["n_genes_min"], 300.0)
        self.assertEqual(th["n_genes_min_derived"], base["n_genes_min"])
        self.assertNotEqual(th["n_genes_min"], th["n_genes_min_derived"])
        # Pass masks are driven by the effective bound: a 290-gene cell fails.
        obs2 = pd.DataFrame({
            "total_counts": [2000.0, 2000.0],
            "n_genes_by_counts": [290.0, 1500.0],
            "pct_counts_mt": [1.0, 1.0],
            "pct_counts_ribo": [10.0, 10.0],
        })
        masks = qct.rna_pass_masks(obs2, th, pct_ribo_max=50.0)
        self.assertFalse(bool(masks["n_genes"][0]))
        self.assertTrue(bool(masks["n_genes"][1]))

    def test_derived_keys_equal_effective_when_no_override(self):
        obs = _rna_obs()
        th = qct.rna_thresholds(obs, **_RNA_KW)
        for k in ("total_counts_min", "total_counts_max", "n_genes_min",
                  "n_genes_max", "pct_counts_mt_max"):
            self.assertEqual(th[k], th[f"{k}_derived"])

    def test_override_below_floor_wins(self):
        obs = _rna_obs()
        th = qct.rna_thresholds(obs, **_RNA_KW, n_genes_min_override=100)
        self.assertEqual(th["n_genes_min"], 100.0)
        # Derived still respects the floor (>= 250).
        self.assertGreaterEqual(th["n_genes_min_derived"], 250.0)


class AtacOverrideMathTests(unittest.TestCase):
    def test_override_sets_applied_and_keeps_derived(self):
        rng = np.random.default_rng(0)
        n_frag = rng.integers(1500, 20000, 500).astype(float)
        f_lo, f_hi, mad_lo, (lo_d, hi_d) = qct.atac_n_fragment_bounds(
            n_frag, k_mad=5.0, n_frag_floor=1500.0, n_fragments_min_override=5000,
        )
        self.assertEqual(f_lo, 5000.0)
        self.assertNotEqual(lo_d, 5000.0)
        self.assertGreaterEqual(lo_d, 1500.0)

    def test_no_override_applied_equals_derived(self):
        rng = np.random.default_rng(1)
        n_frag = rng.integers(1500, 20000, 500).astype(float)
        f_lo, f_hi, mad_lo, (lo_d, hi_d) = qct.atac_n_fragment_bounds(
            n_frag, k_mad=5.0, n_frag_floor=1500.0,
        )
        self.assertEqual(f_lo, lo_d)
        self.assertEqual(f_hi, hi_d)


class OverrideFigureMarkerTests(unittest.TestCase):
    def test_override_is_red_and_mad_is_grey(self):
        # Override 300 displaces the MAD/floor-derived 250.
        markers, _, _ = figures.build_mad_range_markers(
            applied_lo=300.0,
            applied_hi=8000.0,
            default_lo=250.0,
            default_hi=8000.0,
            default_mad_lo_raw=180.0,
            default_floor=250.0,
            hi_skip_above=1_000_000,
            log_axis=True,
            derived_lo=250.0,
        )
        active = {m[0]: m[2] for m in markers}
        # Chosen override drawn red (active), MAD-derived + raw drawn grey.
        self.assertIn(300.0, active)
        self.assertTrue(active[300.0])
        self.assertIn(250.0, active)
        self.assertFalse(active[250.0])
        self.assertIn(180.0, active)
        self.assertFalse(active[180.0])

    def test_no_override_path_unchanged(self):
        kw = dict(
            applied_lo=250.0, applied_hi=8000.0, default_lo=250.0, default_hi=8000.0,
            default_mad_lo_raw=180.0, default_floor=250.0,
            hi_skip_above=1_000_000, log_axis=True,
        )
        before, lo, hi = figures.build_mad_range_markers(**kw)
        after, lo2, hi2 = figures.build_mad_range_markers(**kw, derived_lo=None, derived_hi=None)
        self.assertEqual(before, after)
        self.assertEqual((lo, hi), (lo2, hi2))

    def test_upper_only_override_grey_mad(self):
        markers, _, _ = figures.build_upper_only_markers(
            applied_hi=15.0,
            default_hi=8.0,
            default_mad_hi_raw=8.0,
            hi_skip_above=100.0,
            pct=True,
            default_fixed_refs=[(5.0, "5%"), (10.0, "10%"), (20.0, "20%")],
            derived_hi=8.0,
        )
        active = {round(m[0], 1): m[2] for m in markers}
        self.assertTrue(active[15.0])     # override red
        self.assertFalse(active[8.0])     # MAD-derived grey


class OverrideProvenanceTests(unittest.TestCase):
    def test_user_override_recorded_without_method(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "parameters.yaml"
            # Mirror the stage's conditional recording for an active override.
            _prov.set_param(
                path, "s1_rna_qc.n_genes_min", 300.0,
                source="user", confidence="high",
                rationale="Manual override (was MAD-derived 250)",
            )
            params = _prov.load(path)
            entry = params["s1_rna_qc.n_genes_min"]
            self.assertEqual(entry["source"], "user")
            self.assertIn("Manual override", entry["rationale"])
            self.assertNotIn("method", entry)

    def test_derived_recording_still_requires_method(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "parameters.yaml"
            with self.assertRaises(ValueError):
                _prov.set_param(
                    path, "s1_rna_qc.n_genes_min", 250.0,
                    source="derived", confidence="high", rationale="MAD",
                )


class DoubletThresholdLineTests(unittest.TestCase):
    """save_figure closes the figure, so capture the axes' lines at save time."""

    def _dashed_xs(self, threshold):
        from unittest import mock
        from executor import figures as _fig
        from executor.stages import post_qc_review as pqr
        scores = np.array([0.1, 0.2, 0.3, 0.6, 0.9])
        flags = np.array([False, False, False, True, True])
        captured: list[float] = []
        real_save = _fig.save_figure

        def _capture(fig, *a, **k):
            ax = fig.axes[0]
            captured.extend(
                float(ln.get_xdata()[0]) for ln in ax.get_lines()
                if ln.get_linestyle() == "--"
            )
            return real_save(fig, *a, **k)

        with tempfile.TemporaryDirectory() as d, \
                mock.patch.object(_fig, "save_figure", _capture):
            pqr._plot_score_hist(scores, flags, title="t", out_dir=Path(d),
                                 stem="dub", threshold=threshold)
        return captured

    def test_threshold_draws_one_red_line(self):
        xs = self._dashed_xs(0.3)
        self.assertEqual(len(xs), 1)
        self.assertAlmostEqual(xs[0], 0.3)

    def test_no_threshold_draws_no_cutoff(self):
        self.assertEqual(self._dashed_xs(None), [])


if __name__ == "__main__":
    unittest.main()
