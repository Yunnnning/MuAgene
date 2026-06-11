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
    "n_fragments_min": 1500, "n_fragments_max": 100000,
    "tss_enrichment_min": 1.5, "tss_enrichment_max": 50,
    "nucleosome_signal_max": 3, "frip_min": 0.25,
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
        self.assertIn("| pct_counts_mt | ≤ 20 | 3 |", t)
        self.assertIn("| pct_counts_ribo | ≤ 50 | 1 |", t)
        self.assertIn("| total_counts | 500 – 50000 | 10 |", t)
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
            frip_removed=ATAC_RM["frip_min"],
        )
        self.assertIn("| parameter | value | cells removed |", t)
        self.assertIn("| nucleosome_signal | < 3 | 1 |", t)
        self.assertIn("| frip | ≥ 0.25 | 6 |", t)
        self.assertIn("| n_fragments | 1500 – 100000 | 4 |", t)

    def test_custom_frip_display_runtime(self):
        t = qc_tables.atac_removal_table(
            ATAC_TH, ATAC_RM, value_label="threshold", include_note=True,
            frip_threshold_display="≥ 0.25 _(computed at runtime)_", frip_removed="—",
        )
        self.assertIn("| parameter | threshold | cells removed | note |", t)
        self.assertIn("≥ 0.25 _(computed at runtime)_", t)
        self.assertIn("| frip | ≥ 0.25 _(computed at runtime)_ | — |", t)

    def test_frip_no_peaks_note(self):
        t = qc_tables.atac_removal_table(
            ATAC_TH, ATAC_RM, value_label="value", include_note=False,
            frip_threshold_display="≥ 0.25 _(not applied — no peaks available)_",
            frip_removed="",
        )
        self.assertIn("not applied — no peaks available", t)


if __name__ == "__main__":
    unittest.main()
