"""A plan-review `revise` must reach the stage that runs, not just the rendered plan.

Regression for stages that read the frozen preprocessing_plan.json directly (S1a
method/max_contamination, S4, S5, S6, S8) and for the s1a_ambient_propose rule, which
looked the override up under a nested key parameters.yaml never contains.
"""
from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path

import yaml

from executor import provenance

_ROOT = Path(__file__).resolve().parents[1]


def _params(tmp: str, entries: dict) -> Path:
    p = Path(tmp) / "parameters.yaml"
    p.write_text(yaml.safe_dump(entries))
    return p


class PlannedValueTests(unittest.TestCase):
    PLAN = {"hvg_n_top_genes": {"value": 2000}}

    def test_user_revise_wins(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = _params(tmp, {"s4_rna_norm.hvg_n_top_genes": {"value": 3000, "source": "user"}})
            self.assertEqual(
                provenance.planned_value(p, self.PLAN, "s4_rna_norm", "hvg_n_top_genes"), 3000)

    def test_echoed_derived_value_is_not_fed_back(self):
        # S4 echoes min(cap, 0.1*n_genes) as source=derived; a rerun must start from the cap.
        with tempfile.TemporaryDirectory() as tmp:
            p = _params(tmp, {"s4_rna_norm.hvg_n_top_genes": {"value": 30, "source": "derived"}})
            self.assertEqual(
                provenance.planned_value(p, self.PLAN, "s4_rna_norm", "hvg_n_top_genes"), 2000)

    def test_missing_falls_back_to_plan_then_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "absent.yaml"
            self.assertEqual(provenance.planned_value(p, self.PLAN, "s", "hvg_n_top_genes"), 2000)
            self.assertEqual(provenance.planned_value(p, {}, "s", "x", "d"), "d")


class S1aMethodRevisedToNoneTests(unittest.TestCase):
    """`revise s1a_ambient method=none` at plan review must disable correction."""

    def test_user_none_passes_through(self):
        import anndata as ad
        import numpy as np
        from executor.stages import s1a_ambient

        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp)
            s0 = run / "internal" / "artifacts" / "s0_ingest"
            s0.mkdir(parents=True)
            ad.AnnData(X=np.ones((5, 3), dtype="float32")).write_h5ad(s0 / "rna_ingest.h5ad")
            provenance.set_param(run / "internal" / "parameters.yaml", "s1a_ambient.method",
                                 "none", source="user", confidence="high", rationale="test")
            plan = {"stages": {"s1a_ambient": {"parameters": {"method": {"value": "auto"}}}}}
            out = s1a_ambient.run(run, plan)
            self.assertEqual(out["method"], "none")


class ContextReportDocxTests(unittest.TestCase):
    """The context gate advertises .docx/.pdf reports; they must actually parse."""

    def test_docx_report_is_read_through_doc_ingest(self):
        try:
            import docx
        except ImportError:
            self.skipTest("python-docx unavailable")
        from executor import context

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ctx.docx"
            d = docx.Document()
            for line in ("Organism: Homo sapiens", "Tissue: PBMC", "Assay: 10x Multiome"):
                d.add_paragraph(line)
            d.save(str(path))
            self.assertFalse(context.is_unfilled_template(path))
            parsed = context.parse_report(context.read_report_text(path))
            self.assertEqual(parsed["organism"], "Homo sapiens")


class SnakemakeRulesUseOverlayTests(unittest.TestCase):
    def test_s1a_propose_reads_flat_key(self):
        text = (_ROOT / "workflow" / "rules" / "s1a_ambient.smk").read_text()
        self.assertNotIn('pdata.get("s1a_ambient", {})', text)
        self.assertIn("planned_value", text)

    def test_s1a_execute_is_gated_on_plan_review(self):
        text = (_ROOT / "workflow" / "rules" / "s1a_ambient.smk").read_text()
        execute = text.split("rule s1a_ambient_execute:", 1)[1].split("output:", 1)[0]
        self.assertRegex(execute, re.escape('"plan_review.approved"'))

    def test_downstream_stages_do_not_read_frozen_values_directly(self):
        stages = _ROOT / "executor" / "stages"
        for name in ("s4_rna_norm", "s5_atac_spectral", "s8_umap"):
            text = (stages / f"{name}.py").read_text()
            self.assertNotRegex(text, r'p\["\w+"\]\["value"\]', name)
        s6 = (stages / "s6_neighbors.py").read_text()
        self.assertNotRegex(s6, r's6_params\.get\("\w+", \{\}\)\.get\("value"')
