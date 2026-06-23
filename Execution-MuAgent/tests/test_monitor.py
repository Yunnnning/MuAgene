"""Unit tests for HPC monitor helpers (no scratch dirs in the package tree)."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import inspect

import yaml

from unittest import mock

from execution_muagent.monitor import (
    DEFAULT_CHECK_INTERVAL_S,
    JobHealth,
    MonitorState,
    REPOLL_BUFFER_S,
    SiteConfig,
    StageSpec,
    Submission,
    _is_run_scoped_progress_path,
    _looks_like_hdf5,
    _looks_like_parquet,
    _persist_snapshot,
    _workflow_complete,
    classify_investigation,
    collect_snapshot,
    discover_child_job_ids,
    evaluate_snapshot_definitive,
    extract_error_context,
    monitor_once,
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

    def test_slurm_job_name_includes_run_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "runs" / "whelanC57A"
            repo_root = Path(tmp) / "repo"
            script = self._render(run_dir, repo_root)
            self.assertIn("#SBATCH --job-name=pma_head_job_whelanC57A", script)


    def test_run_name_exported_to_child_jobs(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "runs" / "whelanC57A"
            repo_root = Path(tmp) / "repo"
            script = self._render(run_dir, repo_root)
            self.assertIn("export PMA_RUN_NAME=whelanC57A", script)


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

    def test_check_interval_constant_is_270(self):
        # The CLI --interval defaults reference this constant; keep it 270.0.
        self.assertEqual(DEFAULT_CHECK_INTERVAL_S, 270.0)

    def test_monitor_state_tolerance_default(self):
        self.assertEqual(MonitorState().tolerance_n, 20)


class RepollCadenceTests(unittest.TestCase):
    def test_collect_snapshot_records_repoll_cadence(self):
        # The snapshot carries interval_s + next_recheck_after_s so Processing-MuAgent
        # derives its re-poll wakeup delay instead of hardcoding a magic number.
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            (run_dir / "internal" / "artifacts").mkdir(parents=True)
            sub = _make_submission(run_dir)
            with mock.patch(
                "execution_muagent.monitor.query_scheduler",
                return_value={"scheduler": "slurm", "state": "RUNNING"},
            ):
                snap = collect_snapshot(sub, interval_s=DEFAULT_CHECK_INTERVAL_S)
        self.assertEqual(snap["interval_s"], DEFAULT_CHECK_INTERVAL_S)
        self.assertEqual(
            snap["next_recheck_after_s"], DEFAULT_CHECK_INTERVAL_S + REPOLL_BUFFER_S
        )
        self.assertEqual(snap["next_recheck_after_s"], 295.0)


class HeadJobVerificationTests(unittest.TestCase):
    def test_head_job_with_populated_outputs_is_verified(self):
        # When Processing populates the head_job spec's outputs (planning S0 submit),
        # the monitor verifies them like any stage and emits stage_output_verified.
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            art = run_dir / "internal" / "artifacts" / "s0_ingest"
            art.mkdir(parents=True)
            good = art / "validation_report.json"
            good.write_text('{"ok": true}', encoding="utf-8")
            meta = run_dir / "internal" / "stage_meta"
            meta.mkdir(parents=True, exist_ok=True)
            (meta / "head_job.yaml").write_text(yaml.safe_dump({
                "schema_version": "1",
                "stage": "head_job",
                "science_description": "orchestrator",
                "resources": {"cpus": 1, "mem_mb": 1000, "walltime_min": 60},
                "inputs": {},
                "outputs": {"validation_report": str(good)},
                "progress_timeout_hint": 120,
            }), encoding="utf-8")
            results = verify_stage_outputs(_make_submission(run_dir))
        self.assertIn("head_job", results)
        self.assertTrue(results["head_job"]["validation_report"][0])

    def test_head_job_without_outputs_is_skipped(self):
        # Multi-stage targets leave outputs empty → still skipped (no spurious finding).
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            meta = run_dir / "internal" / "stage_meta"
            meta.mkdir(parents=True, exist_ok=True)
            (meta / "head_job.yaml").write_text(yaml.safe_dump({
                "schema_version": "1",
                "stage": "head_job",
                "science_description": "orchestrator",
                "resources": {"cpus": 1, "mem_mb": 1000, "walltime_min": 60},
                "inputs": {},
                "outputs": {},
                "progress_timeout_hint": 120,
            }), encoding="utf-8")
            results = verify_stage_outputs(_make_submission(run_dir))
        self.assertNotIn("head_job", results)


class ErrorContextTests(unittest.TestCase):
    """The monitor must surface the real child-log exception, not just the generic
    Snakemake envelope, so Processing-MuAgent can root-cause from `hpc-status`."""

    def _run_with_child_log(self, tmp: str, child_text: str, head_text: str):
        run_dir = Path(tmp) / "run"
        slurm_logs = (
            run_dir / "internal" / "snakemake" / ".snakemake" / "slurm_logs"
            / "rule_s5_atac_spectral_execute"
        )
        slurm_logs.mkdir(parents=True)
        (slurm_logs / "1015981.log").write_text(child_text, encoding="utf-8")
        head = run_dir / "head.out"
        head.write_text(head_text, encoding="utf-8")
        sub = Submission(
            agent="Processing-MuAgent", executor="slurm", job_id="1015843",
            run_dir=str(run_dir), config=str(run_dir / "cfg.yaml"),
            target="s7_clustering_propose", repo_root=str(run_dir / "repo"),
            log_path=str(head), submitted_at="2026-06-12T13:46:00Z",
        )
        return sub

    def test_extract_error_context_prefers_child_traceback(self) -> None:
        child = (
            "RuleException:\n"
            'UnboundLocalError in file ".../s5_atac_spectral.smk", line 38:\n'
            "cannot access local variable '_io' where it is not associated with a value\n"
            "2026-06-12 15:09:56 - ERROR - RuleException:\n"
            "WorkflowError:\n"
        )
        head = "Error in rule s5_atac_spectral_execute:\nWorkflowError:\n"
        with tempfile.TemporaryDirectory() as tmp:
            sub = self._run_with_child_log(tmp, child, head)
            ctx = extract_error_context(sub)
        # Names the real exception + the file:line — not just the generic envelope.
        self.assertIn("UnboundLocalError", ctx)
        self.assertIn("s5_atac_spectral.smk", ctx)
        self.assertIn("1015981.log", ctx)

    def test_workflow_error_marker_finding_carries_root_cause(self) -> None:
        child = (
            "RuleException:\n"
            'UnboundLocalError in file ".../s5_atac_spectral.smk", line 38:\n'
        )
        head = "Error in rule s5_atac_spectral_execute:\nWorkflowError:\n"
        with tempfile.TemporaryDirectory() as tmp:
            sub = self._run_with_child_log(tmp, child, head)
            snapshot = collect_snapshot(sub)
            findings = evaluate_snapshot_definitive(snapshot)
        marker = next(f for f in findings if f.code == "workflow_error_marker")
        self.assertIn("Root cause", marker.message)
        self.assertIn("UnboundLocalError", marker.message)

    def test_no_error_context_when_no_markers(self) -> None:
        head = "rule s4_rna_norm_execute:\nFinished jobid: 2\n"
        with tempfile.TemporaryDirectory() as tmp:
            sub = self._run_with_child_log(tmp, "all good\nFinished jobid: 9\n", head)
            snapshot = collect_snapshot(sub)
        self.assertEqual(snapshot.get("error_context"), "")


class KillActionOwnershipTests(unittest.TestCase):
    """kill_action is owned by the in-memory monitor loop, never inherited from a
    prior snapshot file — so a resubmit (fresh MonitorState) cannot show a phantom
    cancelled step from a dead run."""

    def _snapshot_path(self, run_dir: Path) -> Path:
        return run_dir / "internal" / "hpc_monitor" / "latest_snapshot.json"

    def test_fresh_loop_does_not_inherit_prior_kill_action(self) -> None:
        import json as _json
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            sub = _make_submission(run_dir)
            path = self._snapshot_path(run_dir)
            path.parent.mkdir(parents=True, exist_ok=True)
            # Simulate a dead prior run's snapshot carrying a kill_action.
            path.write_text(_json.dumps({
                "submission": {"job_id": "OLD-1015843"},
                "kill_action": {"attempted": True, "head_job_id": "OLD-1015843"},
            }), encoding="utf-8")
            # A fresh loop (killed nothing) persists a finding-free snapshot.
            _persist_snapshot(sub, {"checked_at": "now"}, MonitorState(), findings=[], cancel_result=None)
            written = _json.loads(path.read_text())
            self.assertIsNone(written["kill_action"])

    def test_kill_action_persists_within_loop_via_state(self) -> None:
        import json as _json
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            sub = _make_submission(run_dir)
            path = self._snapshot_path(run_dir)
            kill = {"attempted": True, "head_job_id": sub.job_id}
            # First: the check that kills records kill_action via cancel_result.
            _persist_snapshot(sub, {"checked_at": "t1"},
                              MonitorState(kill_action=kill), findings=[], cancel_result=kill)
            self.assertEqual(_json.loads(path.read_text())["kill_action"], kill)
            # Next finding-free check (cancel_result=None) keeps it via threaded state.
            _persist_snapshot(sub, {"checked_at": "t2"},
                              MonitorState(kill_action=kill), findings=[], cancel_result=None)
            self.assertEqual(_json.loads(path.read_text())["kill_action"], kill)


class DeltaLivenessTests(unittest.TestCase):
    """Liveness is judged by CPU *activity* between investigations, not mere presence —
    so a finished-but-lingering / deadlocked process (flat CPU, still holding memory) is
    correctly classified confirmed_dead instead of 'recovered'. Memory is diagnostic
    only (MaxRSS is monotonic) and not part of the classification."""

    BASE = {
        "scheduler_state": "RUNNING",
        "error_markers": [],
        "child_storage_hang_ids": [],
        "filesystem_responsive": True,
    }

    def _classify(self, **kw):
        return classify_investigation({**self.BASE, **kw})[0]

    def test_flat_cpu_delta_confirmed_dead(self) -> None:
        # Memory held + CPU consumed historically, but CPU did NOT advance → dead.
        self.assertEqual(self._classify(cpu_delta=0.0), "confirmed_dead")

    def test_cpu_advancing_recovered(self) -> None:
        self.assertEqual(self._classify(cpu_delta=42.0), "recovered")

    def test_first_sample_no_prior_recovers_as_baseline(self) -> None:
        # cpu_delta None (no previous investigation) → can't assert activity → recover.
        self.assertEqual(self._classify(cpu_delta=None), "recovered")

    def test_error_marker_still_dead_regardless_of_deltas(self) -> None:
        v, _ = classify_investigation({**self.BASE, "error_markers": ["WorkflowError"],
                                       "cpu_delta": 99.0})
        self.assertEqual(v, "confirmed_dead")


class WorkflowCompleteDetectionTests(unittest.TestCase):
    def test_complete_marker_no_errors_is_true(self) -> None:
        snap = {"error_markers": [], "head_log_tail": "3 of 3 steps (100%) done\nComplete log(s): /x\n",
                "latest_snakemake_tail": ""}
        self.assertTrue(_workflow_complete(snap))

    def test_error_marker_blocks_completion(self) -> None:
        snap = {"error_markers": ["WorkflowError"], "head_log_tail": "3 of 3 steps (100%) done",
                "latest_snakemake_tail": ""}
        self.assertFalse(_workflow_complete(snap))

    def test_partial_progress_is_not_complete(self) -> None:
        snap = {"error_markers": [], "head_log_tail": "2 of 3 steps (67%) done",
                "latest_snakemake_tail": ""}
        self.assertFalse(_workflow_complete(snap))


class MonitorOnceWorkflowCompleteTests(unittest.TestCase):
    """A head job whose workflow finished cleanly but is still RUNNING (lingering)
    is cancelled (head only) and marked DONE — without recording a kill_action."""

    def _synthetic_snapshot(self) -> dict:
        return {
            "scheduler": {"state": "RUNNING"},
            "error_markers": [],
            "head_log_tail": "3 of 3 steps (100%) done\nComplete log(s): /x\n",
            "latest_snakemake_tail": "",
            "latest_progress_file": {"exists": True, "mtime": 1.0},
            "head_log": {"size": 10},
            "child_job_ids": [],
        }

    def test_lingering_head_cancelled_head_only_and_done(self) -> None:
        from execution_muagent import monitor as m
        cancel_calls: list = []
        with tempfile.TemporaryDirectory() as tmp:
            sub = _make_submission(Path(tmp) / "run")
            with mock.patch.object(m, "collect_snapshot", return_value=self._synthetic_snapshot()), \
                 mock.patch.object(m, "verify_stage_outputs", return_value={}), \
                 mock.patch.object(m, "evaluate_snapshot_definitive", return_value=[]), \
                 mock.patch.object(m, "write_report", return_value=None), \
                 mock.patch.object(m, "cancel_submission_jobs",
                                   side_effect=lambda s, ids, t=5: cancel_calls.append(list(ids)) or {}):
                _, findings, cancel_result, _, state = m.monitor_once(
                    sub, MonitorState(), kill_on_hang=True)
        self.assertEqual(state.health, JobHealth.DONE)
        self.assertTrue(any(f.code == "workflow_complete" for f in findings))
        self.assertIsNone(cancel_result)        # not a kill → monitor_once returns no cancel_result
        self.assertIsNone(state.kill_action)     # so Processing sees no phantom cancellation
        self.assertEqual(cancel_calls, [[]])     # head only (empty child-id list)

    def test_kill_on_hang_false_marks_done_without_cancel(self) -> None:
        from execution_muagent import monitor as m
        with tempfile.TemporaryDirectory() as tmp:
            sub = _make_submission(Path(tmp) / "run")
            with mock.patch.object(m, "collect_snapshot", return_value=self._synthetic_snapshot()), \
                 mock.patch.object(m, "verify_stage_outputs", return_value={}), \
                 mock.patch.object(m, "evaluate_snapshot_definitive", return_value=[]), \
                 mock.patch.object(m, "write_report", return_value=None), \
                 mock.patch.object(m, "cancel_submission_jobs") as cancel:
                _, _, _, _, state = m.monitor_once(sub, MonitorState(), kill_on_hang=False)
        self.assertEqual(state.health, JobHealth.DONE)
        cancel.assert_not_called()


class MonitorWatchTerminalBreakTests(unittest.TestCase):
    def test_breaks_on_done_even_if_jobs_still_active(self) -> None:
        from execution_muagent import monitor as m
        done = MonitorState(health=JobHealth.DONE)
        with tempfile.TemporaryDirectory() as tmp:
            sub = _make_submission(Path(tmp) / "run")
            with mock.patch.object(m, "monitor_once", return_value=({}, [], None, None, done)), \
                 mock.patch.object(m, "submission_jobs_active", return_value=True) as sja, \
                 mock.patch.object(m, "validate_terminal_outputs", return_value=[]):
                m.monitor_watch(sub, interval_s=0.01, max_checks=3)
        sja.assert_not_called()  # DONE break precedes the jobs-active re-poll


if __name__ == "__main__":
    unittest.main()
