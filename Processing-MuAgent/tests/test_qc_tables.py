"""Tests for the shared RNA/ATAC QC removal-table builders."""
from __future__ import annotations

import unittest

from executor import qc_tables


RNA_TH = {
    "total_counts_min": 500, "total_counts_max": 50000,
    "n_genes_min": 200, "n_genes_max": 8000,
    "pct_counts_mt_max": 20, "pct_counts_ribo_max": 50,
}
RNA_RM = {
    "total_counts": 10, "n_genes": 5, "pct_counts_mt": 3,
    "pct_counts_ribo": 1, "multiple_metrics": 2, "total_removed": 15,
}

ATAC_TH = {
    "n_fragments_min": 1000, "n_fragments_max": 100000,
    "tss_enrichment_min": 1.5, "tss_enrichment_max": 50,
    "nucleosome_signal_max": 3, "frip_min": 0.2,
}
ATAC_RM = {
    "n_fragments": 4, "tss_enrichment": 2, "nucleosome_signal": 1,
    "frip_min": 6, "multiple_metrics": 1, "total_removed": 12,
}


class RnaTableTests(unittest.TestCase):
    def test_value_label_no_note(self):
        t = qc_tables.rna_removal_table(RNA_TH, RNA_RM, value_label="value", include_note=False)
        self.assertIn("| parameter | value | cells removed |", t)
        self.assertNotIn("note", t.splitlines()[0])
        # Removal-condition format: cells outside [lo, hi] are removed
        self.assertIn("| pct_counts_mt | > 20 | 3 |", t)
        self.assertIn("| pct_counts_ribo | > 50 | 1 |", t)
        self.assertIn("| total_counts | < 500 or > 50000 | 10 |", t)
        self.assertIn("| multiple_metrics | — | 2 |", t)
        self.assertIn("| total_removed | — | 15 |", t)

    def test_threshold_label_with_note(self):
        t = qc_tables.rna_removal_table(RNA_TH, RNA_RM, value_label="threshold", include_note=True)
        self.assertIn("| parameter | threshold | cells removed | note |", t)
        # Note column carries the per-metric explanation.
        self.assertIn(qc_tables.NOTES["pct_counts_mt"], t)


class AtacTableTests(unittest.TestCase):
    def test_default_frip_display(self):
        t = qc_tables.atac_removal_table(
            ATAC_TH, ATAC_RM, value_label="value", include_note=False,
            peak_source="user_peaks",
            frip_removed=ATAC_RM["frip_min"],
        )
        self.assertIn("| parameter | value | cells removed |", t)
        # Removal-condition format
        self.assertIn("| nucleosome_signal | ≥ 3 | 1 |", t)
        self.assertIn("| frip | < 0.20 | 6 |", t)
        self.assertIn("| n_fragments | < 1000 or > 100000 | 4 |", t)

    def test_custom_frip_display_runtime(self):
        t = qc_tables.atac_removal_table(
            ATAC_TH, ATAC_RM, value_label="threshold", include_note=True,
            frip_threshold_display="< 0.2 _(computed at runtime)_", frip_removed="—",
        )
        self.assertIn("| parameter | threshold | cells removed | note |", t)
        self.assertIn("< 0.2 _(computed at runtime)_", t)
        self.assertIn("| frip | < 0.2 _(computed at runtime)_ | — |", t)

    def test_frip_no_peaks_note(self):
        t = qc_tables.atac_removal_table(
            ATAC_TH, ATAC_RM, value_label="value", include_note=False,
            frip_threshold_display="< 0.2 _(not applied — no peaks available)_",
            frip_removed="",
        )
        self.assertIn("not applied — no peaks available", t)


class SkipDisplayTests(unittest.TestCase):
    """Skip-sentinel display helpers produce clean labels for workaround parameters."""

    def test_rna_full_skip_total_counts(self):
        th = {**RNA_TH, "total_counts_min": 0, "total_counts_max": 1e10}
        t = qc_tables.rna_removal_table(th, RNA_RM)
        self.assertIn("| total_counts | not applied |", t)

    def test_rna_upper_skip_total_counts_keeps_floor(self):
        th = {**RNA_TH, "total_counts_min": 500, "total_counts_max": 1e10}
        t = qc_tables.rna_removal_table(th, RNA_RM)
        self.assertIn("| total_counts | < 500 |", t)

    def test_rna_full_skip_n_genes(self):
        th = {**RNA_TH, "n_genes_min": 0, "n_genes_max": 5e6}
        t = qc_tables.rna_removal_table(th, RNA_RM)
        self.assertIn("| n_genes | not applied |", t)

    def test_rna_pct_mt_disabled(self):
        th = {**RNA_TH, "pct_counts_mt_max": 100}
        t = qc_tables.rna_removal_table(th, RNA_RM)
        self.assertIn("| pct_counts_mt | not applied |", t)

    def test_rna_pct_ribo_disabled(self):
        th = {**RNA_TH, "pct_counts_ribo_max": 100}
        t = qc_tables.rna_removal_table(th, RNA_RM)
        self.assertIn("| pct_counts_ribo | not applied |", t)

    def test_atac_n_fragments_full_skip(self):
        th = {**ATAC_TH, "n_fragments_min": 0, "n_fragments_max": 1e9}
        t = qc_tables.atac_removal_table(th, ATAC_RM, frip_removed=0)
        self.assertIn("| n_fragments | not applied |", t)

    def test_atac_tss_fully_disabled(self):
        th = {**ATAC_TH, "tss_enrichment_min": 0, "tss_enrichment_max": 999}
        t = qc_tables.atac_removal_table(th, ATAC_RM, frip_removed=0)
        self.assertIn("| tss_enrichment | not applied |", t)

    def test_atac_tss_upper_only(self):
        th = {**ATAC_TH, "tss_enrichment_min": 1.5, "tss_enrichment_max": 999}
        t = qc_tables.atac_removal_table(th, ATAC_RM, frip_removed=0)
        self.assertIn("| tss_enrichment | < 1.50 |", t)

    def test_atac_nucleosome_disabled(self):
        th = {**ATAC_TH, "nucleosome_signal_max": 999}
        t = qc_tables.atac_removal_table(th, ATAC_RM, frip_removed=0)
        self.assertIn("| nucleosome_signal | not applied |", t)

    def test_atac_frip_disabled(self):
        th = {**ATAC_TH, "frip_min": 0}
        t = qc_tables.atac_removal_table(th, ATAC_RM, frip_removed=0)
        self.assertIn("| frip | not applied |", t)

    def test_atac_frip_disabled_overrides_runtime_note(self):
        th = {**ATAC_TH, "frip_min": 0}
        t = qc_tables.atac_removal_table(
            th, ATAC_RM, frip_runtime_note=True, frip_removed=0,
        )
        self.assertIn("| frip | not applied |", t)
        self.assertNotIn("computed at runtime", t)

    def test_atac_frip_no_peaks_when_enabled(self):
        th = {**ATAC_TH, "frip_min": 0.2}
        t = qc_tables.atac_removal_table(th, ATAC_RM, peak_source=None, frip_removed=0)
        self.assertIn("| frip | < 0.20 _(not applied — no peaks available)_ |", t)

    def test_normal_values_unchanged(self):
        t = qc_tables.rna_removal_table(RNA_TH, RNA_RM)
        self.assertIn("| total_counts | < 500 or > 50000 |", t)
        self.assertIn("| pct_counts_mt | > 20 |", t)


class FlowStepTests(unittest.TestCase):
    def test_preprocessing_flow_steps_split_rna_atac_qc(self):
        from pathlib import Path
        from executor import qc_summary as qcs

        counts = {
            "rna_raw": 100, "atac_raw_barcodes": 100,
            "rna_qc_post": 80, "atac_qc_post": 70,
            "rna_post_doublet": 60, "atac_post_doublet": 60,
            "n_cells_joint": 60,
        }
        steps = qcs._preprocessing_flow_steps(
            Path("/tmp/run"), counts, "paired", include_final_stage=False,
        )
        self.assertEqual([s["stage"] for s in steps], [
            "1. raw",
            "2. after RNA QC",
            "3. after ATAC QC",
            "4. after doublet removal",
        ])
        self.assertEqual(steps[1]["rna"], 80)
        self.assertEqual(steps[1]["atac"], 100)
        self.assertEqual(steps[2]["rna"], 80)
        self.assertEqual(steps[2]["atac"], 70)

    def test_flow_chart_label(self):
        from executor import qc_summary as qcs

        self.assertEqual(qcs._flow_chart_label("2. after RNA QC"), "after RNA\nQC")
        self.assertEqual(qcs._flow_chart_label("1. raw"), "raw")


if __name__ == "__main__":
    unittest.main()
