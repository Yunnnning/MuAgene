"""Regression guard: plan_review.approved must gate every sN_execute rule.

Snakemake rule files are pure Python modulo the rule bodies; we read them as
text and assert the string `plan_review.approved` appears in the input block of
each sN_execute rule. This is a lightweight check — a full Snakemake DAG parse
would also work but takes ~2s and requires snakemake in the test env.
"""
from __future__ import annotations

import re
from pathlib import Path

RULES_DIR = Path(__file__).parent.parent / "workflow" / "rules"

STAGES = ["s1_rna_qc", "s2_atac_qc", "s3_doublets", "s4_rna_norm",
          "s5_atac_lsi", "s6_dimred", "s7_clustering", "s8_umap"]


def _execute_rule_body(stage: str) -> str:
    """Return the text of `rule <stage>_execute:` up to the next `rule` / EOF."""
    text = (RULES_DIR / f"{stage}.smk").read_text()
    m = re.search(rf"rule {stage}_execute:\s*\n(.*?)(?=\nrule |\Z)", text, re.DOTALL)
    assert m, f"Could not locate rule {stage}_execute in {stage}.smk"
    return m.group(1)


def test_plan_review_gates_every_execute_rule() -> None:
    """The plan addendum requires plan_review.approved as an input on every sN_execute."""
    missing: list[str] = []
    for stage in STAGES:
        body = _execute_rule_body(stage)
        if "plan_review.approved" not in body:
            missing.append(stage)
    assert not missing, (
        f"Stages missing `plan_review.approved` input: {missing}. "
        "This is the execution-boundary contract — S1..S8 execute must block on "
        "plan_review approval."
    )
