"""Unit tests for HPC monitor helpers (no scratch dirs in the package tree)."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import inspect

import yaml

from execution_muagent.monitor import (
    MonitorState,
    SiteConfig,
    StageSpec,
    Submission,
    _is_run_scoped_progress_path,
    _looks_like_hdf5,
    _looks_like_parquet,
    discover_child_job_ids,
    monitor_watch,
    parse_job_ids_from_log,
    render_submission_script,
    validate_terminal_outputs,
    verify_output_file,
    verify_stage_outputs,
)

_HDF5_SIG = b"\x89HDF\r\n\x1a\n"


def _make_submission(run_dir: Path) -> Submission:
    return Submission(
        agent="Processing-MuAgent",
        executor="slurm",
        job_id="992574",
        run_dir=str(run_dir),
        config=str(run_dir / "cfg.yaml"),
        target="all",
        repo_root=str(run_dir / "repo"),
        log_path=str(run_dir / "head.out"),
        submitted_at="2026-06-01T11:57:46Z",
    )


class MonitorHelperTests(unittest.TestCase):
    def test_parse_job_ids_snakemake_external_jobid(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "snakemake.log"
            log.write_text(
                "Submitted job 10 with external jobid '992575'.\n"
                "Submitted job 2 with external jobid '992578'.\n",
                encoding="utf-8",
            )
            self.assertEqual(parse_job_ids_from_log(log), ["992575", "992578"])

    def test_parse_job_ids_from_slurm_log_filename(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "992579.log"
            log.write_text("Finished jobid: 0\n", encoding="utf-8")
            self.assertEqual(parse_job_ids_from_log(log), ["992579"])

    def test_discover_child_job_ids_from_run_logs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_dir = root / "run"
            slurm_logs = (
                run_dir / "internal" / "snakemake" / ".snakemake" / "slurm_logs" / "rule_x"
            )
            slurm_logs.mkdir(parents=True)
            (slurm_logs / "992578.log").write_text(
                "Submitted job 2 with external jobid '992578'.\n",
                encoding="utf-8",
            )
            submission = Submission(
                agent="Processing-MuAgent",
                executor="slurm",
                job_id="992574",
                run_dir=str(run_dir),
                config=str(run_dir / "cfg.yaml"),
                target="all",
                repo_root=str(root / "repo"),
                log_path=str(root / "repo" / "logs" / "pma_runner-992574.out"),
                submitted_at="2026-06-01T11:57:46Z",
            )
            self.assertEqual(discover_child_job_ids(submission), ["992578"])

    def test_run_scoped_progress_excludes_other_run_head_log(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_a = root / "run_a"
            run_b = root / "run_b"
            run_a.mkdir()
            run_b.mkdir()
            other_log = root / "logs" / "pma_runner-999.out"
            other_log.parent.mkdir()
            other_log.write_text("other run\n", encoding="utf-8")
            submission = Submission(
                agent="Processing-MuAgent",
                executor="slurm",
                job_id="992574",
                run_dir=str(run_a),
                config=str(run_a / "cfg.yaml"),
                target="all",
                repo_root=str(root / "repo"),
                log_path=str(root / "logs" / "pma_runner-992574.out"),
                submitted_at="2026-06-01T11:57:46Z",
            )
            self.assertFalse(_is_run_scoped_progress_path(submission, other_log))
            self.assertTrue(
                _is_run_scoped_progress_path(submission, run_a / "internal" / "log.jsonl")
            )


def _make_spec(**overrides) -> StageSpec:
    defaults = dict(
        schema_version="1",
        stage="head_job",
        science_description="Orchestrator",
        resources={"cpus": 1, "mem_mb": 4000, "walltime_min": 1440},
        inputs={},
        outputs={},
        progress_timeout_hint=120.0,
    )
    defaults.update(overrides)
    return StageSpec(**defaults)


def _make_slurm_site(*, partition="cpu", account="vaquerizas", conda_env="grn") -> SiteConfig:
    return SiteConfig(
        scheduler="slurm",
        partition=partition,
        account=account,
        conda_env=conda_env,
    )


class RenderSubmissionScriptTests(unittest.TestCase):
    """Regression tests for render_submission_script.

    Bug fixed: PMA_CONFIG was built from repo_root instead of run_dir, and
    launch_runner.sh was called without --configfile / --profile / target args,
    causing `KeyError: 'run_dir'` in the Snakefile on first snakemake invocation.
    """

    def _render(self, run_dir, repo_root, site_config=None, target="s1a_ambient_execute"):
        spec = _make_spec()
        sc = site_config or _make_slurm_site()
        log = Path(run_dir) / "internal" / "hpc_monitor" / "logs" / "head.out"
        return render_submission_script(spec, sc, repo_root, run_dir, log, target)

    def test_pma_config_uses_run_dir_not_repo_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "runs" / "MyRun"
            repo_root = Path(tmp) / "MuAgene" / "Processing-MuAgent"
            script = self._render(run_dir, repo_root)
            # run_dir must appear in PMA_CONFIG
            self.assertIn(str(run_dir.resolve()), script)
            # repo_root must NOT appear in PMA_CONFIG line
            for line in script.splitlines():
                if "PMA_CONFIG=" in line:
                    self.assertNotIn(str(repo_root), line,
                                     "PMA_CONFIG must point at run_dir, not repo_root")

    def test_launch_runner_receives_configfile_arg(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            repo_root = Path(tmp) / "repo"
            script = self._render(run_dir, repo_root)
            self.assertIn("--configfile $PMA_CONFIG", script)

    def test_launch_runner_receives_profile_arg(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            repo_root = Path(tmp) / "repo"
            script = self._render(run_dir, repo_root)
            self.assertIn("--profile", script)
            self.assertIn("profiles/slurm", script)

    def test_launch_runner_receives_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            repo_root = Path(tmp) / "repo"
            script = self._render(run_dir, repo_root, target="s3_doublets_execute")
            self.assertIn("s3_doublets_execute", script)

    def test_slurm_partition_and_account_exported(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            repo_root = Path(tmp) / "repo"
            sc = _make_slurm_site(partition="hmem", account="mylab")
            script = self._render(run_dir, repo_root, site_config=sc)
            self.assertIn("export PMA_SLURM_PARTITION=hmem", script)
            self.assertIn("export PMA_SLURM_ACCOUNT=mylab", script)

    def test_conda_env_exported(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            repo_root = Path(tmp) / "repo"
            sc = _make_slurm_site(conda_env="myenv")
            script = self._render(run_dir, repo_root, site_config=sc)
            self.assertIn("export PMA_CONDA_ENV=myenv", script)

    def test_pbs_queue_and_project_exported(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            repo_root = Path(tmp) / "repo"
            sc = SiteConfig(
                scheduler="pbs", queue="workq", project="vaquerizas", conda_env="grn"
            )
            script = self._render(run_dir, repo_root, site_config=sc)
            self.assertIn("export PMA_PBS_QUEUE=workq", script)
            self.assertIn("export PMA_PBS_PROJECT=vaquerizas", script)
            self.assertIn("profiles/pbs", script)

    def test_run_dir_distinct_from_repo_root(self):
        """Canonical contract: run_dir and repo_root are independent paths."""
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "runs" / "MyExperiment"
            repo_root = Path(tmp) / "code" / "Processing-MuAgent"
            script = self._render(run_dir, repo_root)
            config_line = next(l for l in script.splitlines() if "PMA_CONFIG=" in l)
            # run_dir segment must appear; repo_root code segment must NOT
            self.assertIn("MyExperiment", config_line)
            self.assertNotIn("Processing-MuAgent", config_line)


class OutputVerificationTests(unittest.TestCase):
    """verify_output_file: proper integrity, not just non-empty."""

    def test_missing_and_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing = Path(tmp) / "nope.h5ad"
            self.assertEqual(verify_output_file(missing), (False, "missing"))
            empty = Path(tmp) / "empty.h5ad"
            empty.touch()
            self.assertEqual(verify_output_file(empty), (False, "empty"))

    def test_json_valid_and_corrupt(self):
        with tempfile.TemporaryDirectory() as tmp:
            good = Path(tmp) / "ok.json"
            good.write_text('{"a": 1}', encoding="utf-8")
            self.assertTrue(verify_output_file(good)[0])
            bad = Path(tmp) / "bad.json"
            bad.write_text("{not json", encoding="utf-8")
            self.assertFalse(verify_output_file(bad)[0])

    def test_text_sentinel_nonempty(self):
        with tempfile.TemporaryDirectory() as tmp:
            sentinel = Path(tmp) / "s8_done.txt"
            sentinel.write_text("done", encoding="utf-8")
            ok, reason = verify_output_file(sentinel)
            self.assertTrue(ok)
            self.assertEqual(reason, "nonempty")

    def test_garbage_h5ad_rejected(self):
        # Not a valid HDF5 file: rejected whether or not h5py is installed.
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "rna_qc.h5ad"
            bad.write_bytes(b"this is not hdf5" * 100)
            self.assertFalse(verify_output_file(bad)[0])

    def test_looks_like_hdf5_signature(self):
        with tempfile.TemporaryDirectory() as tmp:
            good = Path(tmp) / "a.h5ad"
            good.write_bytes(_HDF5_SIG + b"\x00" * 64)
            self.assertTrue(_looks_like_hdf5(good))
            bad = Path(tmp) / "b.h5ad"
            bad.write_bytes(b"\x00" * 64)
            self.assertFalse(_looks_like_hdf5(bad))

    def test_looks_like_parquet_magic(self):
        with tempfile.TemporaryDirectory() as tmp:
            good = Path(tmp) / "calls.parquet"
            good.write_bytes(b"PAR1" + b"\x00" * 16 + b"PAR1")
            self.assertTrue(_looks_like_parquet(good))
            bad = Path(tmp) / "x.parquet"
            bad.write_bytes(b"nope" + b"\x00" * 16 + b"nope")
            self.assertFalse(_looks_like_parquet(bad))


class TerminalValidationTests(unittest.TestCase):
    def _write_stage_meta(self, run_dir: Path, stage: str, outputs: dict[str, str]) -> None:
        meta_dir = run_dir / "internal" / "stage_meta"
        meta_dir.mkdir(parents=True, exist_ok=True)
        (meta_dir / f"{stage}.yaml").write_text(yaml.safe_dump({
            "schema_version": "1",
            "stage": stage,
            "science_description": stage,
            "resources": {"cpus": 1, "mem_mb": 1000, "walltime_min": 60},
            "inputs": {},
            "outputs": outputs,
            "progress_timeout_hint": 30,
        }), encoding="utf-8")

    def test_corrupt_output_reported(self):
        # S1 now declares qc_summary.json as its monitored output (not rna_qc.h5ad,
        # which is deleted after post_qc_review approval).
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            run_dir.mkdir()
            art = run_dir / "internal" / "artifacts" / "s1_rna_qc"
            art.mkdir(parents=True)
            bad = art / "qc_summary.json"
            bad.write_bytes(b"not-valid-json{{{")
            self._write_stage_meta(run_dir, "s1_rna_qc", {"qc_summary_json": str(bad)})
            findings = validate_terminal_outputs(_make_submission(run_dir))
            self.assertEqual([f.code for f in findings], ["output_missing"])

    def test_valid_outputs_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            run_dir.mkdir()
            art = run_dir / "internal" / "artifacts" / "s5_atac_spectral"
            art.mkdir(parents=True)
            good = art / "spectral_summary.json"
            good.write_text('{"ok": true}', encoding="utf-8")
            self._write_stage_meta(run_dir, "s5_atac_spectral", {"summary": str(good)})
            self.assertEqual(validate_terminal_outputs(_make_submission(run_dir)), [])
            # verify_stage_outputs reports the stage as verified.
            results = verify_stage_outputs(_make_submission(run_dir))
            self.assertTrue(results["s5_atac_spectral"]["summary"][0])


class DefaultsTests(unittest.TestCase):
    def test_watch_interval_default_is_270(self):
        self.assertEqual(
            inspect.signature(monitor_watch).parameters["interval_s"].default, 270.0
        )

    def test_monitor_state_tolerance_default(self):
        self.assertEqual(MonitorState().tolerance_n, 20)


if __name__ == "__main__":
    unittest.main()
