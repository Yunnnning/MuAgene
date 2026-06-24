"""Tests for _cleanup_qc_intermediates — post-approval h5ad removal."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import yaml

from executor.cli import _cleanup_qc_intermediates
from executor.run_paths import RunPaths


def _init_run(tmp: str, *, run_config: dict | None = None) -> RunPaths:
    paths = RunPaths(tmp)
    paths.ensure()
    paths.parameters_yaml.write_text(
        yaml.safe_dump({"plan": {"workflow_branch": "paired"}})
    )
    if run_config is not None:
        cfg_path = paths.deliv_config / "run.yaml"
        cfg_path.parent.mkdir(parents=True, exist_ok=True)
        cfg_path.write_text(yaml.safe_dump(run_config))
    return paths


def _touch(path: Path, content: bytes = b"placeholder") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


class CleanupQCIntermediatesTests(unittest.TestCase):
    def test_deletes_target_h5ads(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            rna_qc   = _touch(paths.artifact("s1_rna_qc",  "rna_qc.h5ad"))
            atac_qc  = _touch(paths.artifact("s2_atac_qc", "atac_qc.h5ad"))
            atac_snap = _touch(paths.artifact("s2_atac_qc", "atac_snap.h5ad"))
            atac_snap_explore = _touch(paths.artifact("qc_explore", "atac_snap_explore.h5ad"))

            deleted = _cleanup_qc_intermediates(Path(tmp))

            self.assertFalse(rna_qc.exists())
            self.assertFalse(atac_qc.exists())
            self.assertFalse(atac_snap.exists())
            self.assertFalse(atac_snap_explore.exists())
            self.assertEqual(len(deleted), 4)

    def test_returns_only_paths_that_existed(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            # Only create two of the three targets
            _touch(paths.artifact("s1_rna_qc",  "rna_qc.h5ad"))
            _touch(paths.artifact("s2_atac_qc", "atac_qc.h5ad"))
            # atac_snap.h5ad is absent (already cleaned or never created)

            deleted = _cleanup_qc_intermediates(Path(tmp))
            self.assertEqual(len(deleted), 2)

    def test_preserves_qc_summary_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            _touch(paths.artifact("s1_rna_qc",  "rna_qc.h5ad"))
            _touch(paths.artifact("s2_atac_qc", "atac_qc.h5ad"))
            _touch(paths.artifact("s2_atac_qc", "atac_snap.h5ad"))
            s1_json  = _touch(paths.artifact("s1_rna_qc",  "qc_summary.json"), b"{}")
            s2_json  = _touch(paths.artifact("s2_atac_qc", "qc_summary.json"), b"{}")

            _cleanup_qc_intermediates(Path(tmp))

            self.assertTrue(s1_json.exists())
            self.assertTrue(s2_json.exists())

    def test_preserves_qc_metrics_parquets(self):
        """The s1/s2 qc_metrics parquets must survive cleanup: they are the durable
        QC-metrics record consumed by the QC-review summary."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            _touch(paths.artifact("s1_rna_qc", "rna_qc.h5ad"))
            kept = [
                _touch(paths.artifact("s1_rna_qc", "qc_metrics_pre.parquet"),  b"PAR1\x00PAR1"),
                _touch(paths.artifact("s1_rna_qc", "qc_metrics_post.parquet"), b"PAR1\x00PAR1"),
                _touch(paths.artifact("s2_atac_qc", "qc_metrics_pre.parquet"),  b"PAR1\x00PAR1"),
                _touch(paths.artifact("s2_atac_qc", "qc_metrics_post.parquet"), b"PAR1\x00PAR1"),
            ]

            _cleanup_qc_intermediates(Path(tmp))

            for p in kept:
                self.assertTrue(p.exists(), f"Expected {p} to be preserved")

    def _make_fragment_caches(self, paths: RunPaths) -> list[Path]:
        caches: list[Path] = []
        for stage in ("qc_explore", "s2_atac_qc"):
            for name in ("atac_fragments_cbf_chrnorm.tsv.gz", "atac_fragments_cbf.tsv.gz"):
                caches.append(_touch(paths.artifact(stage, name), b"data"))
                caches.append(_touch(paths.artifact(stage, name + ".tbi"), b"idx"))
        return caches

    def test_retains_cbf_fragment_caches_by_default(self):
        """Default (retain_for_integration unset/true): the chr-normalised fragment
        caches are KEPT past the gate — Integration-MuAgent re-counts them against a
        consensus peak set, reading their contents, not just the recorded filename."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)  # no run.yaml -> default retain
            caches = self._make_fragment_caches(paths)

            deleted = _cleanup_qc_intermediates(Path(tmp))

            for p in caches:
                self.assertTrue(p.exists(), f"Expected {p} to be retained")
                self.assertNotIn(str(p), deleted)

    def test_retains_cbf_fragment_caches_when_flag_true(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, run_config={"retain_for_integration": True})
            caches = self._make_fragment_caches(paths)

            _cleanup_qc_intermediates(Path(tmp))

            for p in caches:
                self.assertTrue(p.exists(), f"Expected {p} to be retained")

    def test_deletes_cbf_fragment_caches_when_opted_out(self):
        """A single-sample-and-done run can set retain_for_integration: false to
        delete the (large) fragment caches and reclaim disk."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp, run_config={"retain_for_integration": False})
            caches = self._make_fragment_caches(paths)

            deleted = _cleanup_qc_intermediates(Path(tmp))

            for p in caches:
                self.assertFalse(p.exists(), f"Expected {p} to be deleted")
            self.assertEqual(len(deleted), len(caches))

    def test_deletes_s1a_recompute_caches(self):
        """S1a recompute caches are dead after approval (no S1a re-run can occur)."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            caches = [
                _touch(paths.artifact("s1a_ambient", "tsne_coords_cache.parquet"), b"PAR1"),
                _touch(paths.artifact("s1a_ambient", "cell_totals.parquet"), b"PAR1"),
            ]

            _cleanup_qc_intermediates(Path(tmp))

            for p in caches:
                self.assertFalse(p.exists(), f"Expected {p} to be deleted")

    def test_preserves_s1a_provenance_diagnostics(self):
        """Per-cell ambient provenance is kept (not a cache, negligible size)."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            kept = [
                _touch(paths.artifact("s1a_ambient", "contamination.parquet"), b"PAR1"),
                _touch(paths.artifact("s1a_ambient", "marker_gene_check.json"), b"{}"),
                _touch(paths.artifact("s1a_ambient", "summary.json"), b"{}"),
            ]

            _cleanup_qc_intermediates(Path(tmp))

            for p in kept:
                self.assertTrue(p.exists(), f"Expected {p} to be preserved")

    def test_preserves_s3_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            kept = [
                _touch(paths.artifact("s3_doublets", "rna_post_doublet.h5ad")),
                _touch(paths.artifact("s3_doublets", "atac_post_doublet.h5ad")),
                _touch(paths.artifact("s3_doublets", "calls.parquet"), b"PAR1\x00PAR1"),
            ]

            _cleanup_qc_intermediates(Path(tmp))

            for p in kept:
                self.assertTrue(p.exists(), f"Expected {p} to be preserved")

    def test_deletes_s0_s1a_heavy_rna_caches(self):
        """rna_ingest.h5ad (~200 MB), metadata_minimal.tsv, and rna_decontaminated.h5ad
        (~400 MB) are content-dead once QC is approved (S4+ read the post-QC h5mu)."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            doomed = [
                _touch(paths.artifact("s0_ingest", "rna_ingest.h5ad")),
                _touch(paths.artifact("s0_ingest", "metadata_minimal.tsv")),
                _touch(paths.artifact("s1a_ambient", "rna_decontaminated.h5ad")),
            ]

            deleted = _cleanup_qc_intermediates(Path(tmp))

            for p in doomed:
                self.assertFalse(p.exists(), f"Expected {p} to be deleted")
            for p in doomed:
                self.assertIn(str(p), deleted)

    def test_preserves_s0_s1a_markers_after_cleanup(self):
        """The durable markers (s0 validation_report.json, s1a summary.json) must
        survive — they carry status + the DAG edges. validation_report.json is also
        read post-gate by S5."""
        from executor.stage_progress import execute_done
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            _touch(paths.artifact("s0_ingest", "rna_ingest.h5ad"))
            _touch(paths.artifact("s1a_ambient", "rna_decontaminated.h5ad"))
            report = _touch(paths.artifact("s0_ingest", "validation_report.json"), b"{}")
            summary = _touch(paths.artifact("s1a_ambient", "summary.json"), b"{}")

            _cleanup_qc_intermediates(Path(tmp))

            self.assertTrue(report.exists())
            self.assertTrue(summary.exists())
            # s1a stays done off summary.json (its EXECUTE_MARKER) after cleanup.
            self.assertTrue(execute_done(paths, "s1a_ambient"))

    def test_preserves_qc_explore_metric_parquets(self):
        """The per-cell QC metric parquets (under qc_explore/) must survive cleanup
        so a post-approval `revise` can re-derive thresholds without a heavy reload."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = _init_run(tmp)
            _touch(paths.artifact("s1_rna_qc", "rna_qc.h5ad"))
            _touch(paths.artifact("s2_atac_qc", "atac_qc.h5ad"))
            kept = [
                _touch(paths.artifact("qc_explore", "rna_qc_metrics.parquet"), b"PAR1\x00PAR1"),
                _touch(paths.artifact("qc_explore", "atac_qc_metrics.parquet"), b"PAR1\x00PAR1"),
                _touch(paths.artifact("qc_explore", "qc_explore.json"), b"{}"),
            ]

            _cleanup_qc_intermediates(Path(tmp))

            for p in kept:
                self.assertTrue(p.exists(), f"Expected {p} to be preserved")

    def test_no_targets_present_returns_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            _init_run(tmp)
            deleted = _cleanup_qc_intermediates(Path(tmp))
            self.assertEqual(deleted, [])


class QcCleanupCommandTests(unittest.TestCase):
    """The standalone `executor qc-cleanup` command: gated on post_qc_review approval,
    deletes the same set as the approve-time cleanup."""

    def _cfg(self, tmp: str) -> str:
        paths = _init_run(tmp)
        cfg = paths.deliv_config / "run.yaml"
        cfg.parent.mkdir(parents=True, exist_ok=True)
        cfg.write_text(yaml.safe_dump({"run_dir": tmp}))
        return str(cfg)

    def test_refuses_when_not_approved(self):
        from click.testing import CliRunner
        from executor.cli import main
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self._cfg(tmp)
            target = _touch(RunPaths(tmp).artifact("s1a_ambient", "rna_decontaminated.h5ad"))
            res = CliRunner().invoke(main, ["qc-cleanup", "--config", cfg])
            self.assertNotEqual(res.exit_code, 0)
            self.assertIn("not approved", res.output)
            self.assertTrue(target.exists(), "must not delete while QC unapproved")

    def test_cleans_when_approved(self):
        from click.testing import CliRunner
        from executor.cli import main
        with tempfile.TemporaryDirectory() as tmp:
            cfg = self._cfg(tmp)
            paths = RunPaths(tmp)
            paths.approved_sentinel("post_qc_review").write_text("")
            doomed = [
                _touch(paths.artifact("s0_ingest", "rna_ingest.h5ad")),
                _touch(paths.artifact("s1a_ambient", "rna_decontaminated.h5ad")),
            ]
            kept = _touch(paths.artifact("s1a_ambient", "summary.json"), b"{}")

            res = CliRunner().invoke(main, ["qc-cleanup", "--config", cfg])

            self.assertEqual(res.exit_code, 0, res.output)
            for p in doomed:
                self.assertFalse(p.exists())
            self.assertTrue(kept.exists())


class QCMatrixUntrackedTests(unittest.TestCase):
    """The QC matrices (rna_qc.h5ad / atac_qc.h5ad) must NOT be declared Snakemake
    rule outputs, and S3 must depend on the durable qc_summary.json marker instead.
    This is the structural fix that lets _cleanup_qc_intermediates delete the
    matrices without making any rule's declared output "missing" — so a
    post-approval `submit` does not re-run S1/S2/S3. (Replaces the earlier temp()
    approach.)
    """

    _RULES = Path(__file__).resolve().parents[1] / "workflow" / "rules"

    @staticmethod
    def _code(text: str) -> str:
        """Drop comment lines so assertions test declarations, not prose."""
        return "\n".join(l for l in text.splitlines() if not l.strip().startswith("#"))

    def _output_block(self, rule_file: str, rule_name: str) -> str:
        """Return the code (comments stripped) of `rule_name`'s `output:` block."""
        import re
        text = (self._RULES / rule_file).read_text()
        sub = text[text.index(f"rule {rule_name}:"):]
        m = re.search(r"\n    output:\n(.*?)\n    (?:params|threads|resources):",
                      sub, re.DOTALL)
        return self._code(m.group(1) if m else "")

    def test_s1_rna_qc_h5ad_not_declared_output(self):
        block = self._output_block("s1_rna_qc.smk", "s1_rna_qc_execute")
        self.assertNotIn("rna_qc.h5ad", block,
                         "s1_rna_qc must NOT declare rna_qc.h5ad as an output")
        self.assertIn("qc_summary.json", block,
                      "s1_rna_qc must declare qc_summary.json as its output")

    def test_s2_atac_qc_h5ad_not_declared_output(self):
        block = self._output_block("s2_atac_qc.smk", "s2_atac_qc_execute")
        self.assertNotIn("atac_qc.h5ad", block,
                         "s2_atac_qc must NOT declare atac_qc.h5ad as an output")
        self.assertIn("qc_summary.json", block,
                      "s2_atac_qc must declare qc_summary.json as its output")

    def test_no_temp_wrappers_remain(self):
        for rule in ("s1_rna_qc.smk", "s2_atac_qc.smk"):
            text = self._code((self._RULES / rule).read_text())
            self.assertNotIn("temp(", text, f"{rule}: temp() approach is superseded")

    def test_s3_depends_on_qc_summary_not_h5ad(self):
        code = self._code((self._RULES / "s3_doublets.smk").read_text())
        # The s3 input function must reference the durable markers, not the matrices.
        self.assertIn("s1_rna_qc\" / \"qc_summary.json", code)
        self.assertIn("s2_atac_qc\" / \"qc_summary.json", code)
        self.assertNotIn("rna_qc.h5ad", code,
                         "s3 must not declare rna_qc.h5ad as an input edge")
        self.assertNotIn("atac_qc.h5ad", code,
                         "s3 must not declare atac_qc.h5ad as an input edge")

    def test_cleanup_targets_are_not_declared_outputs(self):
        """Every file deleted at the gate must be untracked by Snakemake — if any
        were a declared output, deleting it would make a rule's output "missing"
        and trigger a re-run on the next submit."""
        names = [
            "atac_fragments_cbf_chrnorm.tsv.gz", "atac_fragments_cbf.tsv.gz",
            "atac_snap_explore.h5ad", "atac_snap.h5ad",
            "tsne_coords_cache.parquet", "cell_totals.parquet",
            # S0/S1a heavy RNA caches now deleted at the gate — must be untracked.
            "rna_ingest.h5ad", "metadata_minimal.tsv", "rna_decontaminated.h5ad",
        ]
        joined = "\n".join(self._code(p.read_text()) for p in self._RULES.glob("*.smk"))
        for name in names:
            self.assertNotIn(name, joined,
                             f"{name} must not be a declared Snakemake output")

    def test_s0_s1a_edges_use_durable_markers(self):
        """S1a must depend on s0's validation_report.json (not rna_ingest.h5ad), and
        S1 on s1a's summary.json (not rna_decontaminated.h5ad) — so deleting those
        heavy caches at the gate never makes a declared output 'missing'."""
        s1a = self._code((self._RULES / "s1a_ambient.smk").read_text())
        self.assertIn("validation_report.json", s1a)
        self.assertIn("summary.json", s1a)        # s1a's own declared output marker
        self.assertNotIn("rna_ingest.h5ad", s1a)
        self.assertNotIn("rna_decontaminated.h5ad", s1a)

        s1 = self._code((self._RULES / "s1_rna_qc.smk").read_text())
        self.assertIn("s1a_ambient\" / \"summary.json", s1)
        self.assertNotIn("rna_decontaminated.h5ad", s1)

        s0 = self._code((self._RULES / "s0_ingest.smk").read_text())
        self.assertNotIn("rna_ingest.h5ad", s0,
                         "s0 must not declare rna_ingest.h5ad as an output")


class S2TempSweepTests(unittest.TestCase):
    """_sweep_stage_temps removes leaked S2 scratch without touching real outputs."""

    def test_sweeps_only_scratch(self):
        try:
            from executor.stages.s2_atac_qc import _sweep_stage_temps
        except ModuleNotFoundError as e:  # numpy/snapatac2 absent in a bare env
            self.skipTest(f"s2_atac_qc dependencies unavailable: {e}")
        with tempfile.TemporaryDirectory() as tmp:
            art = Path(tmp)
            scratch = [
                _touch(art / "tmpABCDEF.h5ad"),
                _touch(art / "_frip_tmp.h5ad"),
                _touch(art / "_peaks_stripped_tmp.bed"),
            ]
            kept = [
                _touch(art / "atac_qc.h5ad"),
                _touch(art / "qc_summary.json", b"{}"),
                _touch(art / "qc_metrics_post.parquet", b"PAR1"),
            ]

            _sweep_stage_temps(art)

            for p in scratch:
                self.assertFalse(p.exists(), f"Expected scratch {p} to be swept")
            for p in kept:
                self.assertTrue(p.exists(), f"Expected {p} to be preserved")


if __name__ == "__main__":
    unittest.main()
