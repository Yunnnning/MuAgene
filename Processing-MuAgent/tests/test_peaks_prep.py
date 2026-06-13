"""Peak-BED preparation for SnapATAC2 (executor/io.prepare_peaks_for_snapatac).

Covers the S5 crash where a Cell Ranger ARC `atac_peaks.bed` (50 `#` comment
lines + Ensembl chrom names) was passed raw to `make_peak_matrix`, which panicked
in Rust (`MissingStartPosition`, a BaseException) and — past the panic — would have
mismatched the UCSC-renamed fragments. The shared helper strips comments and applies
the SAME Ensembl→UCSC convention as the fragment path so peaks and fragments align.
"""
import tempfile
import unittest
from pathlib import Path

from executor import io


class PreparePeaksTests(unittest.TestCase):
    def _prep(self, text: str, **kw) -> list[str]:
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "peaks.bed"
            out = Path(tmp) / "out.bed"
            src.write_text(text)
            io.prepare_peaks_for_snapatac(src, out, **kw)
            return out.read_text().splitlines()

    def test_strips_comments_and_renames_ensembl_to_ucsc(self) -> None:
        text = (
            "# id=Mouse1\n"
            "# pipeline_version=cellranger-arc-2.0.2\n"
            "\n"
            "1\t100\t200\n"
            "X\t10\t20\n"
            "MT\t5\t15\n"
            "GL456211.1\t1\t9\n"     # unplaced scaffold → unchanged
        )
        lines = self._prep(text, add_chr_prefix=True)
        self.assertEqual(lines, [
            "chr1\t100\t200",
            "chrX\t10\t20",
            "chrM\t5\t15",           # MT → chrM (not chrMT)
            "GL456211.1\t1\t9",      # scaffold passed through, matching fragments
        ])

    def test_no_double_prefix_for_already_ucsc_peaks(self) -> None:
        # MACS3 peaks called on renamed fragments are already chr-prefixed.
        lines = self._prep("chr1\t100\t200\nchr2\t5\t9\n", add_chr_prefix=True)
        self.assertEqual(lines, ["chr1\t100\t200", "chr2\t5\t9"])

    def test_no_rename_when_add_chr_prefix_false(self) -> None:
        lines = self._prep("# hdr\n1\t100\t200\n", add_chr_prefix=False)
        self.assertEqual(lines, ["1\t100\t200"])

    def test_handles_space_delimited_and_drops_malformed(self) -> None:
        lines = self._prep("1 100 200\nchrJUNK\nbad\tline\n", add_chr_prefix=True)
        self.assertEqual(lines, ["chr1\t100\t200"])  # space→tab; short/non-numeric dropped

    def test_raises_when_no_valid_intervals(self) -> None:
        with self.assertRaises(ValueError):
            self._prep("# only comments\n\n", add_chr_prefix=True)


if __name__ == "__main__":
    unittest.main()
