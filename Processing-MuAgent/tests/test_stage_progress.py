import tempfile
import unittest
from pathlib import Path

import yaml

from executor import cli
from executor.run_paths import RunPaths
from executor.stage_progress import (
    EXECUTE_MARKERS,
    MONITOR_PIPELINE,
    collect_failed_snakemake_rules,
    collect_monitor_killed_monitor_ids,
    infer_resume_target,
    required_human_approvals,
    snakemake_rules_for_monitor,
    stage_states,
)


def _states_by_label(paths: RunPaths) -> dict[str, str]:
    return {label: state for label, _task, state in stage_states(paths)}


def _write_cluster_rule_log(paths: RunPaths, rule_name: str, text: str) -> None:
    log_dir = paths.snakemake_workdir / ".snakemake" / "slurm_logs" / f"rule_{rule_name}"
    log_dir.mkdir(parents=True, exist_ok=True)
    (log_dir / "1.log").write_text(text)


class StageProgressTests(unittest.TestCase):
    def _mark_planning_done(self, paths: RunPaths) -> None:
        """Seed the merged-planning markers (validation_report + qc_explore JSON)
        so infer_resume_target moves past the planning phase."""
        for stage, marker in (
            ("s0_ingest", "validation_report.json"),
            ("qc_explore", "qc_explore.json"),
        ):
            p = paths.artifact(stage, marker)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("{}")

    def _init_run(self, tmp: str, *, branch: str = "paired",
                  planning_done: bool = True) -> RunPaths:
        paths = RunPaths(tmp)
        paths.ensure()
        paths.parameters_yaml.write_text(
            yaml.safe_dump({"plan": {"workflow_branch": branch}})
        )
        # The merged s0_ingest planning job runs before checkpoint #1; most tests
        # exercise the QC phase and beyond, so treat planning as done by default.
        if planning_done:
            self._mark_planning_done(paths)
        return paths

    def test_every_monitor_step_has_snakemake_rules_except_resolution_gate(self):
        for mid in MONITOR_PIPELINE:
            rules = snakemake_rules_for_monitor(mid)
            if mid == "resolution_review":
                self.assertEqual(rules, ())
            elif mid in EXECUTE_MARKERS or mid in (
                "plan_review", "post_qc_review", "s7_sweep", "s7_labels",
            ):
                self.assertTrue(rules, msg=mid)

    def test_automated_steps_include_propose_and_execute_rules(self):
        self.assertEqual(
            snakemake_rules_for_monitor("s3_doublets"),
            ("s3_doublets_propose", "s3_doublets_execute"),
        )
        self.assertEqual(
            snakemake_rules_for_monitor("s7_sweep"),
            ("s7_clustering_propose",),
        )
        self.assertEqual(
            snakemake_rules_for_monitor("s7_labels"),
            ("s7_clustering_execute",),
        )

    def test_stage_states_lists_substeps_and_resolution_review(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            paths.approved_sentinel("plan_review").write_text("")
            paths.artifact("s1a_ambient", "rna_decontaminated.h5ad").parent.mkdir(
                parents=True, exist_ok=True
            )
            paths.artifact("s1a_ambient", "rna_decontaminated.h5ad").write_text("")
            paths.awaiting_sentinel("post_qc_review").write_text("")

            states = _states_by_label(paths)

            self.assertEqual(states["S1a"], "done")
            self.assertEqual(states["qc_review"], "awaiting_approval")
            self.assertIn("resolution_review", states)

    def test_s1_s2_show_done_after_cleanup(self):
        """qc_summary.json persists after post_qc_review cleanup; S1/S2 must stay done."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            paths.approved_sentinel("plan_review").write_text("")
            paths.approved_sentinel("post_qc_review").write_text("")
            for stage, marker in (
                ("s1a_ambient", "rna_decontaminated.h5ad"),
                ("s1_rna_qc", "qc_summary.json"),
                ("s2_atac_qc", "qc_summary.json"),
                ("s3_doublets", "calls.parquet"),
            ):
                p = paths.artifact(stage, marker)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text("{}")
            states = _states_by_label(paths)
            self.assertEqual(states["S1"], "done")
            self.assertEqual(states["S2"], "done")

    def test_resolution_review_awaiting_after_sweep(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            sweep = paths.artifact("s7_clustering", "sweep.parquet")
            sweep.parent.mkdir(parents=True, exist_ok=True)
            sweep.write_text("")
            paths.awaiting_sentinel("s7_clustering").write_text("")

            states = _states_by_label(paths)
            self.assertEqual(states["S7"], "done")
            self.assertEqual(states["resolution_review"], "awaiting_approval")
            self.assertEqual(states["S7-labels"], "pending")

    def test_s3_shows_failed_from_cluster_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            paths.approved_sentinel("plan_review").write_text("")
            for stage, marker in (
                ("s1a_ambient", "rna_decontaminated.h5ad"),
                ("s1_rna_qc", "qc_summary.json"),
                ("s2_atac_qc", "qc_summary.json"),
            ):
                p = paths.artifact(stage, marker)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text("")

            _write_cluster_rule_log(
                paths, "s3_doublets_execute", "Error in rule s3_doublets_execute:\n",
            )

            self.assertIn("s3_doublets_execute", collect_failed_snakemake_rules(paths))
            self.assertEqual(_states_by_label(paths)["S3"], "failed")

    def test_propose_failure_marks_step_failed(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            paths.approved_sentinel("plan_review").write_text("")
            _write_cluster_rule_log(
                paths, "s2_atac_qc_execute", "Error in rule s2_atac_qc_execute:\n",
            )

            states = _states_by_label(paths)
            self.assertEqual(states["S2"], "failed")

    def test_downstream_blocked_after_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            paths.approved_sentinel("plan_review").write_text("")
            for stage, marker in (
                ("s1a_ambient", "rna_decontaminated.h5ad"),
                ("s1_rna_qc", "qc_summary.json"),
                ("s2_atac_qc", "qc_summary.json"),
            ):
                p = paths.artifact(stage, marker)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text("")
            _write_cluster_rule_log(
                paths, "s3_doublets_execute", "Error in rule s3_doublets_execute:\n",
            )

            states = _states_by_label(paths)
            self.assertEqual(states["S3"], "failed")
            self.assertEqual(states["qc_review"], "blocked")
            self.assertEqual(states["S4"], "blocked")
            self.assertEqual(states["S8"], "blocked")

    def test_done_outputs_override_stale_failure_in_main_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            paths.approved_sentinel("plan_review").write_text("")
            p = paths.artifact("s1a_ambient", "rna_decontaminated.h5ad")
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("")
            _write_cluster_rule_log(
                paths, "s1a_ambient_execute", "Error in rule s1a_ambient_execute:\n",
            )
            smk_log = paths.snakemake_workdir / ".snakemake" / "log"
            smk_log.mkdir(parents=True, exist_ok=True)
            (smk_log / "run.snakemake.log").write_text(
                "Error in rule s1a_ambient_execute:\n"
                "Exiting because a job execution failed.\n"
            )

            self.assertEqual(_states_by_label(paths)["S1a"], "done")

    def test_infer_resume_at_failed_mid_qc_step(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            paths.approved_sentinel("plan_review").write_text("")
            for stage, marker in (
                ("s1a_ambient", "rna_decontaminated.h5ad"),
                ("s1_rna_qc", "qc_summary.json"),
            ):
                p = paths.artifact(stage, marker)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text("")

            # _last_incomplete_execute returns the LAST incomplete stage so Snakemake
            # can chain s2→s3 in one submission; s3 is also missing here.
            self.assertEqual(infer_resume_target(tmp), "s3_doublets_execute")

    def test_infer_resume_after_qc_approval_at_s5(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            paths.approved_sentinel("plan_review").write_text("")
            paths.approved_sentinel("post_qc_review").write_text("")
            for stage, marker in (
                ("s1a_ambient", "rna_decontaminated.h5ad"),
                ("s1_rna_qc", "qc_summary.json"),
                ("s2_atac_qc", "qc_summary.json"),
                ("s3_doublets", "calls.parquet"),
                ("s4_rna_norm", "rna_norm.h5ad"),
            ):
                p = paths.artifact(stage, marker)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text("")

            # _last_incomplete_execute returns the LAST incomplete stage so Snakemake
            # can chain s5→s6 in one submission; s6 is also missing here.
            self.assertEqual(infer_resume_target(tmp), "s6_neighbors_execute")

    def test_infer_resume_targets_planning_first(self):
        """A fresh run (no planning artifacts) routes the merged planning job first,
        so `submit` dispatches it via Execution-MuAgent."""
        with tempfile.TemporaryDirectory() as tmp:
            self._init_run(tmp, planning_done=False)
            self.assertEqual(infer_resume_target(tmp), "s0_ingest_execute")

    def test_infer_resume_plan_review_after_planning(self):
        """Once planning artifacts exist but plan_review is unapproved, the resume
        target is the plan_review render gate."""
        with tempfile.TemporaryDirectory() as tmp:
            self._init_run(tmp, planning_done=True)
            self.assertEqual(infer_resume_target(tmp), "plan_review_propose")

    def test_infer_resume_incomplete_planning_when_explore_missing(self):
        """A planning job that died after the validation report but before the QC
        exploration is still treated as incomplete."""
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp, planning_done=False)
            report = paths.artifact("s0_ingest", "validation_report.json")
            report.parent.mkdir(parents=True, exist_ok=True)
            report.write_text("{}")
            self.assertEqual(infer_resume_target(tmp), "s0_ingest_execute")

    def test_required_human_approvals_planning_targets_empty(self):
        """Planning targets run before checkpoint #1 — they cannot require the
        plan_review sentinel (which would deadlock `submit`)."""
        self.assertEqual(required_human_approvals("s0_ingest_execute"), ())
        self.assertEqual(required_human_approvals("plan_review_propose"), ())
        # QC-phase targets still require plan_review approval.
        self.assertEqual(required_human_approvals("s1a_ambient_execute"), ("plan_review",))

    def test_cli_status_wrapper_matches_stage_progress(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            self.assertEqual(cli._stage_states(paths), stage_states(paths))

    def _write_monitor_kill_report(self, paths: RunPaths, child_job_ids: list[str]) -> None:
        monitor_dir = paths.run_dir / "internal" / "hpc_monitor"
        monitor_dir.mkdir(parents=True, exist_ok=True)
        actions = [
            {"job_id": jid, "role": "child", "attempted": True, "returncode": 0}
            for jid in child_job_ids
        ]
        actions.append({"job_id": "999000", "role": "head", "attempted": True, "returncode": 0})
        payload = {
            "attempted": True,
            "executor": "slurm",
            "head_job_id": "999000",
            "child_job_ids": child_job_ids,
            "actions": actions,
        }
        report = (
            "# HPC Monitor Report\n\n"
            "## Kill Action\n\n"
            f"```json\n{__import__('json').dumps(payload, indent=2)}\n```\n"
        )
        (monitor_dir / "latest_report.md").write_text(report)

    def test_monitor_kill_shows_cancelled_without_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            paths.approved_sentinel("plan_review").write_text("")
            for stage, marker in (
                ("s1a_ambient", "rna_decontaminated.h5ad"),
                ("s1_rna_qc", "qc_summary.json"),
                ("s2_atac_qc", "qc_summary.json"),
            ):
                p = paths.artifact(stage, marker)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text("{}")
            _write_cluster_rule_log(paths, "s3_doublets_execute", "slurmstepd: CANCELLED\n")
            self._write_monitor_kill_report(paths, ["42"])
            (paths.snakemake_workdir / ".snakemake" / "slurm_logs" / "rule_s3_doublets_execute").mkdir(
                parents=True, exist_ok=True
            )
            (paths.snakemake_workdir / ".snakemake" / "slurm_logs" / "rule_s3_doublets_execute" / "42.log").write_text(
                "slurmstepd: CANCELLED\n"
            )

            self.assertIn("s3_doublets", collect_monitor_killed_monitor_ids(paths))
            states = _states_by_label(paths)
            self.assertEqual(states["S3"], "cancelled")
            self.assertEqual(states["qc_review"], "blocked")

    def test_monitor_kill_does_not_override_done_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            paths.approved_sentinel("plan_review").write_text("")
            p = paths.artifact("s1a_ambient", "rna_decontaminated.h5ad")
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("")
            _write_cluster_rule_log(paths, "s1a_ambient_execute", "slurmstepd: CANCELLED\n")
            self._write_monitor_kill_report(paths, ["7"])
            log_dir = paths.snakemake_workdir / ".snakemake" / "slurm_logs" / "rule_s1a_ambient_execute"
            log_dir.mkdir(parents=True, exist_ok=True)
            (log_dir / "7.log").write_text("slurmstepd: CANCELLED\n")

            self.assertEqual(_states_by_label(paths)["S1a"], "done")


if __name__ == "__main__":
    unittest.main()
