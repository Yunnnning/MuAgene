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

    def test_injects_conda_activation_after_shebang(self):
        text = "#!/bin/sh\npython -m snakemake --mode remote\n"
        out = hpc.sanitize_snakemake_jobscript_text(text)
        lines = out.splitlines()
        # Shebang stays first; activation block is injected right after it.
        self.assertEqual(lines[0], "#!/bin/sh")
        self.assertIn(hpc.PMA_ACTIVATION_MARKER, out)
        self.assertIn('conda activate "$PMA_CONDA_ENV"', out)
        # Guarded on PMA_CONDA_ENV so local/non-conda jobs are unaffected.
        self.assertIn('if [ -n "${PMA_CONDA_ENV:-}" ]; then', out)

    def test_conda_activation_injection_is_idempotent(self):
        text = "#!/bin/sh\npython -m snakemake\n"
        once = hpc.sanitize_snakemake_jobscript_text(text)
        twice = hpc.sanitize_snakemake_jobscript_text(once)
        self.assertEqual(once, twice)
        self.assertEqual(once.count(hpc.PMA_ACTIVATION_MARKER), 1)


class ManifestJobIdTests(unittest.TestCase):
    """A resumed submit must report the NEW head-job's id, never a stale entry
    left in execution_manifest.jsonl by a previous submission."""

    def _write_manifest(self, run_dir: Path, job_ids: list[str]) -> None:
        import json
        mp = run_dir / "internal" / "hpc_monitor" / "execution_manifest.jsonl"
        mp.parent.mkdir(parents=True, exist_ok=True)
        mp.write_text("".join(json.dumps({"job_id": j}) + "\n" for j in job_ids))

    def test_entry_count_counts_nonempty_lines(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            self.assertEqual(hpc._manifest_entry_count(run_dir), 0)
            self._write_manifest(run_dir, ["1001", "1002"])
            self.assertEqual(hpc._manifest_entry_count(run_dir), 2)

    def test_wait_ignores_stale_entry_and_times_out(self):
        # Manifest already holds a prior submission's entry; no new entry is
        # appended, so the wait must NOT return that stale id — it must time out.
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            self._write_manifest(run_dir, ["1015641"])
            baseline = hpc._manifest_entry_count(run_dir)
            got = hpc._wait_for_manifest(run_dir, baseline_count=baseline, timeout_s=0.2)
            self.assertIsNone(got)

    def test_wait_returns_new_entry(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            self._write_manifest(run_dir, ["1015641"])
            baseline = hpc._manifest_entry_count(run_dir)
            self._write_manifest(run_dir, ["1015641", "1015843"])
            got = hpc._wait_for_manifest(run_dir, baseline_count=baseline, timeout_s=0.2)
            self.assertEqual(got, "1015843")


if __name__ == "__main__":
    unittest.main()
