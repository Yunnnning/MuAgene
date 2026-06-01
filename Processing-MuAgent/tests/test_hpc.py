from pathlib import Path
import tempfile
import unittest

from executor import hpc


class HpcTests(unittest.TestCase):
    def test_sanitize_snakemake_jobscript_text_removes_storage_local_copies(self):
        text = (
            "python -m snakemake --shared-fs-usage persistence storage-local-copies "
            "input-output --mode remote\n"
        )

        sanitized = hpc.sanitize_snakemake_jobscript_text(text)

        self.assertNotIn("storage-local-copies", sanitized)
        self.assertIn("--shared-fs-usage persistence", sanitized)
        self.assertRegex(sanitized, r"--mode\s+('subprocess'|subprocess)")
        self.assertNotRegex(sanitized, r"--mode\s+('remote'|remote)")

    def test_sanitize_snakemake_jobscript_updates_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            jobscript = Path(tmp) / "snakejob.sh"
            jobscript.write_text(
                "#!/bin/sh\n"
                "python -m snakemake --shared-fs-usage sources storage-local-copies "
                "source-cache --mode remote\n"
            )

            changed = hpc.sanitize_snakemake_jobscript(jobscript)

            self.assertIs(changed, True)
            self.assertNotIn("storage-local-copies", jobscript.read_text())

    def test_sanitize_snakemake_jobscript_rewrites_remote_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            jobscript = Path(tmp) / "snakejob.sh"
            jobscript.write_text(
                "#!/bin/sh\n"
                "python -m snakemake --shared-fs-usage sources source-cache --mode remote\n"
            )

            changed = hpc.sanitize_snakemake_jobscript(jobscript)

            self.assertIs(changed, True)
            self.assertRegex(jobscript.read_text(), r"--mode\s+('subprocess'|subprocess)")


if __name__ == "__main__":
    unittest.main()
