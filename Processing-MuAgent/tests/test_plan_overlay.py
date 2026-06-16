"""Overlay + regeneration: a `revise` must win over the frozen plan everywhere.

The frozen ``preprocessing_plan.json`` is the *default* layer; ``parameters.yaml``
(a user ``revise``) is the *effective* layer. These tests pin the single overlay
rule (``provenance.effective_value`` / ``effective_params`` /
``plan_assembler.overlay_plan``), the S1 fix that previously ignored overrides,
and the ``revise`` behaviour split between the plan_review gate (regenerate the
plan deliverables) and the post_qc_review gate (invalidate downstream).
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import yaml
from click.testing import CliRunner

from executor import plan_assembler, provenance
from executor.cli import main
from executor.run_paths import RunPaths
from executor.stages import s1_rna_qc


def _write_params(path: Path, entries: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(entries))


def _user_entry(value, rationale="user revise"):
    return {"value": value, "source": "user", "confidence": "high",
            "rationale": rationale, "assumptions": []}


class EffectiveValueTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.params_path = Path(self._tmp.name) / "parameters.yaml"

    def tearDown(self):
        self._tmp.cleanup()

    def test_override_wins_over_plan(self):
        _write_params(self.params_path, {"s1_rna_qc.pct_mt_ceiling": _user_entry(15.0)})
        plan_params = {"pct_mt_ceiling": {"value": 20.0}}
        self.assertEqual(
            provenance.effective_value(self.params_path, plan_params, "s1_rna_qc", "pct_mt_ceiling"),
            15.0)

    def test_falls_back_to_plan_then_default(self):
        _write_params(self.params_path, {})
        plan_params = {"pct_mt_ceiling": {"value": 20.0}}
        self.assertEqual(
            provenance.effective_value(self.params_path, plan_params, "s1_rna_qc", "pct_mt_ceiling"),
            20.0)
        self.assertEqual(
            provenance.effective_value(self.params_path, {}, "s1_rna_qc", "missing", 7),
            7)

    def test_effective_params_overlays_and_keeps_others(self):
        _write_params(self.params_path, {"s2_atac_qc.frip_min": _user_entry(0.15)})
        plan_params = {"frip_min": {"value": 0.2}, "tss_enrichment_min": {"value": 1.5}}
        eff = provenance.effective_params(self.params_path, plan_params, "s2_atac_qc")
        self.assertEqual(eff["frip_min"]["value"], 0.15)
        self.assertEqual(eff["tss_enrichment_min"]["value"], 1.5)


class S1HonorsReviseTests(unittest.TestCase):
    """Regression: S1 used to read recipe knobs straight from the frozen plan,
    silently ignoring `revise s1_rna_qc.<knob>`. It must now overlay."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.params_path = Path(self._tmp.name) / "parameters.yaml"

    def tearDown(self):
        self._tmp.cleanup()

    def test_resolve_param_picks_up_override(self):
        _write_params(self.params_path, {"s1_rna_qc.pct_mt_ceiling": _user_entry(12.0)})
        plan_params = {"pct_mt_ceiling": {"value": 20.0}}
        self.assertEqual(
            s1_rna_qc._resolve_param(self.params_path, plan_params, "pct_mt_ceiling", 20.0),
            12.0)

    def test_resolve_param_default_when_absent(self):
        _write_params(self.params_path, {})
        self.assertEqual(
            s1_rna_qc._resolve_param(self.params_path, {}, "pct_mt_ceiling", 20.0),
            20.0)


class OverlayPlanTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.params_path = Path(self._tmp.name) / "parameters.yaml"
        self.plan = {
            "workflow_branch": "paired",
            "stages": {
                "s2_atac_qc": {"parameters": {
                    "frip_min": {"value": 0.2, "source": "default", "rationale": "plan default"},
                    "tss_enrichment_min": {"value": 1.5, "source": "default", "rationale": "keep"},
                }},
            },
        }

    def tearDown(self):
        self._tmp.cleanup()

    def test_user_override_reflected_with_source(self):
        _write_params(self.params_path, {"s2_atac_qc.frip_min": _user_entry(0.15, "loosen FRiP")})
        eff = plan_assembler.overlay_plan(self.plan, self.params_path)
        p = eff["stages"]["s2_atac_qc"]["parameters"]["frip_min"]
        self.assertEqual(p["value"], 0.15)
        self.assertEqual(p["source"], "user")
        self.assertEqual(p["rationale"], "loosen FRiP")
        # Non-overridden parameter untouched.
        self.assertEqual(eff["stages"]["s2_atac_qc"]["parameters"]["tss_enrichment_min"]["value"], 1.5)
        # Original plan not mutated (deep copy).
        self.assertEqual(self.plan["stages"]["s2_atac_qc"]["parameters"]["frip_min"]["value"], 0.2)

    def test_echoed_value_syncs_but_keeps_plan_rationale(self):
        # After a stage runs it echoes the applied value back as source!=user.
        _write_params(self.params_path, {
            "s2_atac_qc.frip_min": {"value": 0.15, "source": "recommended", "confidence": "medium",
                                     "rationale": "stage echo", "assumptions": []}})
        eff = plan_assembler.overlay_plan(self.plan, self.params_path)
        p = eff["stages"]["s2_atac_qc"]["parameters"]["frip_min"]
        self.assertEqual(p["value"], 0.15)            # value synced → display == applied
        self.assertEqual(p["rationale"], "plan default")  # explanatory text preserved


def _seed_run(tmp: str) -> tuple[RunPaths, Path]:
    paths = RunPaths(tmp)
    paths.ensure()
    _write_params(paths.parameters_yaml, {"plan": {"value": "paired", "source": "user",
                                                   "confidence": "high", "rationale": "x"}})
    plan = {
        "workflow_branch": "paired",
        "execution": {"mode": "local", "settings": {}},
        "stages": {
            "s2_atac_qc": {"parameters": {
                "frip_min": {"value": 0.2, "source": "default", "rationale": "plan default"},
            }},
        },
    }
    paths.preprocessing_plan.parent.mkdir(parents=True, exist_ok=True)
    paths.preprocessing_plan.write_text(json.dumps(plan))  # plan is read via json.loads
    qc = paths.artifact("qc_explore", "qc_explore.json")
    qc.parent.mkdir(parents=True, exist_ok=True)
    qc.write_text("{}")
    cfg = Path(tmp) / "run.yaml"
    cfg.write_text(yaml.safe_dump({"run_dir": str(paths.run_dir)}))
    return paths, cfg


class ReviseGateBranchTests(unittest.TestCase):
    """`revise` regenerates plan deliverables before plan approval, and invalidates
    downstream QC artifacts after it."""

    def test_plan_review_unapproved_regenerates_and_keeps_artifacts(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths, cfg = _seed_run(tmp)
            # A QC artifact that must NOT be deleted at the plan-review gate.
            s2_summary = paths.artifact("s2_atac_qc", "qc_summary.json")
            s2_summary.parent.mkdir(parents=True, exist_ok=True)
            s2_summary.write_text("{}")

            res = CliRunner().invoke(main, [
                "revise", "s2_atac_qc", "s2_atac_qc.frip_min=0.15",
                "--config", str(cfg), "--rationale", "loosen FRiP"])
            self.assertEqual(res.exit_code, 0, res.output)
            self.assertIn("Regenerated plan deliverables", res.output)
            # Deliverable regenerated and overlay reflects the override.
            self.assertTrue(paths.plan_review_md.exists())
            self.assertIn("0.15", paths.plan_review_md.read_text())
            # No downstream invalidation at this gate.
            self.assertTrue(s2_summary.exists())
            self.assertNotIn("Invalidated", res.output)

    def test_plan_review_approved_invalidates_downstream(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths, cfg = _seed_run(tmp)
            # Plan approved → post_qc_review gate is the active one.
            approved = paths.approved_sentinel("plan_review")
            approved.parent.mkdir(parents=True, exist_ok=True)
            approved.write_text("ok")
            # Downstream QC artifacts + gate outputs that must be invalidated.
            s2_summary = paths.artifact("s2_atac_qc", "qc_summary.json")
            s2_summary.parent.mkdir(parents=True, exist_ok=True)
            s2_summary.write_text("{}")
            gate = paths.awaiting_sentinel("post_qc_review")
            gate.parent.mkdir(parents=True, exist_ok=True)
            gate.write_text("")

            res = CliRunner().invoke(main, [
                "revise", "s2_atac_qc", "s2_atac_qc.frip_min=0.15",
                "--config", str(cfg), "--rationale", "loosen FRiP"])
            self.assertEqual(res.exit_code, 0, res.output)
            self.assertIn("Invalidated", res.output)
            self.assertFalse(s2_summary.exists())
            self.assertFalse(gate.exists())
            self.assertNotIn("Regenerated plan deliverables", res.output)


if __name__ == "__main__":
    unittest.main()
