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
        metrics = {
            "n_fragments": {"values": rng.lognormal(8, 0.7, 500), "lo": 1500,
                            "hi": 60000, "log": True, "mad_hi": True},
            "tss_enrichment": {"values": rng.gamma(4, 1.2, 500), "lo": 1.5, "hi": 50.0},
            "nucleosome_signal": {"values": rng.gamma(2, 0.4, 500), "lo": None, "hi": 3.0},
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

    def test_pct_mt_shows_mad_ref_when_below_floor(self):
        rng = np.random.default_rng(0)
        vals = rng.uniform(0, 8, 500)
        fig, axes = plt.subplots(1, 1)
        figures.plot_qc_threshold_histograms(
            {
                "pct_counts_mt": {
                    "values": vals,
                    "lo": None,
                    "hi": 5.0,
                    "refs": [(3.4, "3.4% (MAD)"), (5.0, "5%"), (10.0, "10%")],
                },
            },
            out_dir=tempfile.gettempdir(),
            stem="_pct_mt_test",
            title="test",
        )
        # Active hi at 5% + MAD ref at 3.4% + 10% ref; floor ref skipped (== hi).
        fig2, ax = plt.subplots()
        figures._draw_threshold_markers(
            ax,
            [(3.4, "3.4% (MAD)", False), (5.0, "5%", True), (10.0, "10%", False)],
            x_range=20.0,
            pct=True,
        )
        self.assertEqual(len(ax.lines), 3)
        ys = figures._stagger_threshold_label_ys(
            [3.4, 5.0, 10.0], 20.0, active=[False, True, False],
        )
        self.assertEqual(ys[1], figures.THRESHOLD_LABEL_Y_TOP)
        plt.close(fig2)


if __name__ == "__main__":
    unittest.main()
