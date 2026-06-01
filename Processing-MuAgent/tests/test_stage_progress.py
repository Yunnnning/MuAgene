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
    infer_resume_target,
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
    def _init_run(self, tmp: str, *, branch: str = "paired") -> RunPaths:
        paths = RunPaths(tmp)
        paths.ensure()
        paths.parameters_yaml.write_text(
            yaml.safe_dump({"plan": {"workflow_branch": branch}})
        )
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
                ("s1_rna_qc", "rna_qc.h5ad"),
                ("s2_atac_qc", "atac_qc.h5ad"),
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
                ("s1_rna_qc", "rna_qc.h5ad"),
                ("s2_atac_qc", "atac_qc.h5ad"),
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
                ("s1_rna_qc", "rna_qc.h5ad"),
            ):
                p = paths.artifact(stage, marker)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text("")

            self.assertEqual(infer_resume_target(tmp), "s2_atac_qc_execute")

    def test_infer_resume_after_qc_approval_at_s5(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            paths.approved_sentinel("plan_review").write_text("")
            paths.approved_sentinel("post_qc_review").write_text("")
            for stage, marker in (
                ("s1a_ambient", "rna_decontaminated.h5ad"),
                ("s1_rna_qc", "rna_qc.h5ad"),
                ("s2_atac_qc", "atac_qc.h5ad"),
                ("s3_doublets", "calls.parquet"),
                ("s4_rna_norm", "rna_norm.h5ad"),
            ):
                p = paths.artifact(stage, marker)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text("")

            self.assertEqual(infer_resume_target(tmp), "s5_atac_spectral_execute")

    def test_cli_status_wrapper_matches_stage_progress(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = self._init_run(tmp)
            self.assertEqual(cli._stage_states(paths), stage_states(paths))


if __name__ == "__main__":
    unittest.main()
