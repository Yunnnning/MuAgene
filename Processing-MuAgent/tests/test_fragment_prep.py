"""Chromosome-naming reconciliation for ATAC fragments (executor/io.py).

Covers the bug where Cell Ranger ARC Ensembl-named fragments ("1", "MT") never
got renamed to the SnapATAC2 UCSC reference ("chr1", "chrM") because the rename
was coupled to bgzip/tabix being on PATH. Forces those binaries to be "missing"
to exercise the pure-Python gzip fallback and asserts the rename still happens.
"""
import gzip
import tempfile
import unittest
from pathlib import Path

from executor import io


class _FakeGenome:
    """Stand-in for a SnapATAC2 Genome — only `.chrom_sizes` is read."""
    def __init__(self, chrom_sizes):
        self.chrom_sizes = chrom_sizes


# UCSC-named reference (matches SnapATAC2 built-in genomes). No scaffold entry.
UCSC_GENOME = _FakeGenome({"chr1": 1000, "chrM": 200})

# Ensembl-named fragments (no chr prefix), mito "MT", plus an unplaced scaffold.
ENSEMBL_FRAGMENTS = (
    "1\t100\t200\tAAAA\t1\n"        # in bounds  -> chr1
    "1\t950\t1200\tBBBB\t1\n"       # end > size -> dropped
    "MT\t10\t50\tCCCC\t1\n"         # mito       -> chrM
    "MT\t150\t300\tDDDD\t1\n"       # end > size -> dropped
    "GL456211.1\t5\t60\tEEEE\t1\n"  # scaffold   -> passthrough
)


def _write_frags(tmp: Path, text: str) -> Path:
    p = tmp / "fragments.tsv.gz"
    with gzip.open(p, "wt") as f:
        f.write(text)
    return p


def _read_rows(path: Path):
    with gzip.open(path, "rt") as f:
        return [ln.rstrip("\n").split("\t") for ln in f if ln.strip()]


class FragmentPrepTests(unittest.TestCase):
    def test_peek_detects_ensembl_no_prefix(self):
        with tempfile.TemporaryDirectory() as tmp:
            frags = _write_frags(Path(tmp), ENSEMBL_FRAGMENTS)
            _, has_chr = io.peek_fragment_chrom_naming(frags)
            self.assertFalse(has_chr)

    def test_rename_and_bounds_without_bgzip_tabix(self):
        # Force the env-tool-free fallback path regardless of the test machine.
        orig_which = io.shutil.which
        io.shutil.which = lambda name: None if name in ("bgzip", "tabix") else orig_which(name)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                tmpp = Path(tmp)
                frags = _write_frags(tmpp, ENSEMBL_FRAGMENTS)
                out_path, add_chr_prefix = io.prepare_fragments_for_snapatac(
                    frags, UCSC_GENOME, out_dir=tmpp,
                )
                self.assertTrue(add_chr_prefix)
                # chrnorm encoded in the filename so a non-renamed stale file is never reused.
                self.assertTrue(out_path.name.endswith("_cbf_chrnorm.tsv.gz"))
                # No tabix index when tabix is unavailable.
                self.assertFalse(Path(str(out_path) + ".tbi").exists())

                rows = _read_rows(out_path)
                chroms = [r[0] for r in rows]
                self.assertEqual(chroms, ["chr1", "chrM", "GL456211.1"])
                # Out-of-bounds fragments removed; mito mapped MT -> chrM (not chrMT).
                self.assertNotIn("chrMT", chroms)
                self.assertEqual(len(rows), 3)
        finally:
            io.shutil.which = orig_which

    def test_no_rename_when_conventions_match(self):
        # Fragments already UCSC-named + UCSC genome -> add_chr_prefix False.
        with tempfile.TemporaryDirectory() as tmp:
            tmpp = Path(tmp)
            frags = _write_frags(tmpp, "chr1\t100\t200\tAAAA\t1\n")
            out_path, add_chr_prefix = io.prepare_fragments_for_snapatac(
                frags, UCSC_GENOME, out_dir=tmpp,
            )
            self.assertFalse(add_chr_prefix)
            self.assertTrue(out_path.name.endswith("_cbf.tsv.gz"))
            self.assertFalse(out_path.name.endswith("_cbf_chrnorm.tsv.gz"))


if __name__ == "__main__":
    unittest.main()
