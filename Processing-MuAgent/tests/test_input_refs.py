import json
import tempfile
import unittest
from pathlib import Path

from executor import io as _io


class InputRefTests(unittest.TestCase):
    def test_write_and_resolve_symlink(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "raw_feature_bc_matrix.h5"
            source.write_bytes(b"fake-h5")
            ref_base = tmp / "artifacts" / "s0_ingest" / "rna_raw"

            _io.write_input_ref(ref_base, source, fmt="10x_h5")

            self.assertTrue(ref_base.is_symlink())
            self.assertEqual(ref_base.resolve(), source.resolve())
            sidecar = Path(str(ref_base) + ".json")
            self.assertTrue(sidecar.exists())
            meta = json.loads(sidecar.read_text())
            self.assertEqual(meta["format"], "10x_h5")
            self.assertEqual(meta["role"], "rna_raw")

            path, fmt = _io.resolve_input_ref(ref_base)
            self.assertEqual(path, source.resolve())
            self.assertEqual(fmt, "10x_h5")

    def test_has_input_ref_detects_symlink_and_legacy(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            art = tmp / "s0_ingest"
            art.mkdir(parents=True)

            self.assertFalse(_io.has_input_ref(art, "rna_raw"))

            source = tmp / "raw.h5ad"
            source.write_bytes(b"x")
            _io.write_input_ref(art / "rna_raw", source, fmt="h5ad")
            self.assertTrue(_io.has_input_ref(art, "rna_raw"))

            (art / "rna_raw").unlink()
            (art / "rna_raw.json").unlink()
            legacy = art / "rna_raw.h5ad"
            legacy.write_bytes(b"legacy")
            self.assertTrue(_io.has_input_ref(art, "rna_raw"))


if __name__ == "__main__":
    unittest.main()
