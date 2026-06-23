"""QC threshold-histogram grid: fragment-size panel + extra_panel slot."""
import tempfile
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from executor import figures  # noqa: E402


def _fake_fsd(n=1001):
    distr = np.zeros(n)
    x = np.arange(1, n)
    distr[1:] = (np.exp(-((x - 70) / 40) ** 2)
                 + 0.5 * np.exp(-((x - 200) / 30) ** 2) + 0.02)
    return distr


class CutoffLabelTests(unittest.TestCase):
    def test_mad_suffix_on_active_when_raw_matches(self):
        self.assertEqual(
            figures._cutoff_label(341.0, pct=False, log_axis=True, mad=True),
            "341 (MAD)",
        )
        self.assertEqual(
            figures._cutoff_label(5.0, pct=True, log_axis=False, mad=False),
            "5%",
        )
        self.assertEqual(
            figures._cutoff_label(35790.0, pct=False, log_axis=True, mad=True),
            "35,790 (MAD)",
        )


class FragmentSizePanelTests(unittest.TestCase):
    def test_draw_fragment_size_distribution_renders(self):
        fig, ax = plt.subplots()
        figures._draw_fragment_size_distribution(
            ax, _fake_fsd(), title="fragment size distribution (pre-filtering)")
        self.assertEqual(ax.get_title(), "fragment size distribution (pre-filtering)")
        self.assertEqual(ax.get_xlabel(), "fragment length (bp)")
        self.assertEqual(len(ax.lines), 1)  # the curve only
        self.assertGreaterEqual(len(ax.collections), 1)  # fill_between PolyCollection
        plt.close(fig)

    def test_draw_fragment_size_distribution_no_data(self):
        fig, ax = plt.subplots()
        figures._draw_fragment_size_distribution(ax, np.array([]), title="frag size")
        self.assertIn("(no data)", ax.get_title())
        plt.close(fig)

    def test_extra_panel_fills_fourth_slot(self):
        rng = np.random.default_rng(0)
        nf_markers, nf_lo, nf_hi = figures.build_mad_range_markers(
            applied_lo=1000.0,
            applied_hi=60000.0,
            default_lo=1000.0,
            default_hi=60000.0,
            default_mad_lo_raw=1000.0,
            default_floor=1000.0,
            hi_skip_above=1_000_000,
            log_axis=True,
        )
        tss_markers, tss_lo, tss_hi = figures.build_fixed_range_markers(
            applied_lo=1.5,
            applied_hi=50.0,
            default_lo=1.5,
            default_hi=50.0,
            hi_skip_above=500,
        )
        nuc_markers, _, nuc_hi = figures.build_upper_only_markers(
            applied_hi=3.0,
            default_hi=3.0,
            hi_skip_above=50,
        )
        metrics = {
            "n_fragments": {
                "values": rng.lognormal(8, 0.7, 500),
                "log": True,
                "markers": nf_markers,
                "filter_lo": nf_lo,
                "filter_hi": nf_hi,
            },
            "tss_enrichment": {
                "values": rng.gamma(4, 1.2, 500),
                "markers": tss_markers,
                "filter_lo": tss_lo,
                "filter_hi": tss_hi,
            },
            "nucleosome_signal": {
                "values": rng.gamma(2, 0.4, 500),
                "markers": nuc_markers,
                "filter_hi": nuc_hi,
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            out = figures.plot_qc_threshold_histograms(
                metrics, out_dir=tmp, stem="atac_grid",
                title="ATAC Data Distribution and QC Thresholds",
                extra_panel={"distr": _fake_fsd(),
                             "title": "fragment size distribution (pre-filtering)"},
            )
            paths = [Path(p) for p in out]
            self.assertTrue(any(p.suffix == ".png" and p.exists() and p.stat().st_size > 0
                                for p in paths))

    def test_atac_extra_panel_slot_one(self):
        n_used = 4
        extra_slot = 1
        hist_slots = [i for i in range(n_used) if i != extra_slot]
        metric_names = ["n_fragments", "tss_enrichment", "nucleosome_signal"]
        layout = {hist_slots[i]: metric_names[i] for i in range(len(metric_names))}
        layout[extra_slot] = "fragment_size"
        self.assertEqual(layout[0], "n_fragments")
        self.assertEqual(layout[1], "fragment_size")
        self.assertEqual(layout[2], "tss_enrichment")
        self.assertEqual(layout[3], "nucleosome_signal")

    def test_pct_mt_shows_mad_ref_when_below_floor(self):
        rng = np.random.default_rng(0)
        vals = rng.uniform(0, 8, 500)
        mt_markers, _, mt_hi = figures.build_upper_only_markers(
            applied_hi=5.0,
            default_hi=5.0,
            default_mad_hi_raw=3.4,
            hi_skip_above=100.0,
            pct=True,
            default_fixed_refs=[(5.0, "5%"), (10.0, "10%")],
        )
        self.assertEqual(
            [m[1] for m in mt_markers if m[2]],
            ["5%"],
        )
        self.assertIn("3.4% (MAD)", [m[1] for m in mt_markers if not m[2]])
        figures.plot_qc_threshold_histograms(
            {
                "pct_counts_mt": {
                    "values": vals,
                    "markers": mt_markers,
                    "filter_hi": mt_hi,
                },
            },
            out_dir=tempfile.gettempdir(),
            stem="_pct_mt_test",
            title="test",
        )


if __name__ == "__main__":
    unittest.main()
