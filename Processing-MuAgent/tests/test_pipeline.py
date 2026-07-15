"""Unit tests for executor.pipeline — the single source of truth for stage topology.

These lock the per-branch membership against the Snakemake DAG (workflow/rules/*.smk):
in particular that s6_neighbors/s7_clustering/s8_umap run on EVERY branch (no branch
guard in the DAG), which is the drift this module was created to prevent.
"""
from __future__ import annotations

import pytest

from executor import pipeline


def test_branches_are_the_four_known_workflow_branches():
    assert set(pipeline.STAGES_BY_BRANCH) == {"paired", "unpaired", "rna_only", "atac_only"}


def test_paired_and_unpaired_run_the_full_modality_pipeline():
    full = set(pipeline.PIPELINE_STAGE_ORDER)
    assert pipeline.stages_for_branch("paired") == full
    assert pipeline.stages_for_branch("unpaired") == full


def test_legacy_separate_branch_is_rejected():
    with pytest.raises(ValueError, match="Unknown workflow_branch"):
        pipeline.stages_for_branch("separate")


def test_rna_only_drops_atac_stages_keeps_s6_s7():
    rna = pipeline.stages_for_branch("rna_only")
    assert {"s2_atac_qc", "s5_atac_spectral"}.isdisjoint(rna)
    assert {"s1a_ambient", "s1_rna_qc", "s4_rna_norm", "s6_neighbors", "s7_clustering", "s8_umap"} <= rna


def test_atac_only_drops_rna_stages_but_keeps_s6_s7_s8():
    """The DAG runs s6/s7/s8 for atac_only (clustering the ATAC spectral embedding)."""
    atac = pipeline.stages_for_branch("atac_only")
    assert {"s1a_ambient", "s1_rna_qc", "s4_rna_norm"}.isdisjoint(atac)
    assert {"s2_atac_qc", "s5_atac_spectral", "s3_doublets",
            "s6_neighbors", "s7_clustering", "s8_umap"} <= atac


def test_s3_runs_on_every_branch():
    for branch in pipeline.STAGES_BY_BRANCH:
        assert "s3_doublets" in pipeline.stages_for_branch(branch)


def test_stages_for_branch_rejects_unknown():
    with pytest.raises(ValueError):
        pipeline.stages_for_branch("bogus")


def test_human_checkpoints_and_automated_split():
    assert pipeline.HUMAN_CHECKPOINTS == ("plan_review", "post_qc_review")
    assert set(pipeline.AUTOMATED_STAGES).isdisjoint(pipeline.HUMAN_CHECKPOINTS)
    assert set(pipeline.AUTOMATED_STAGES) | set(pipeline.HUMAN_CHECKPOINTS) == set(pipeline.STAGES)


def test_stage_aliasing():
    assert pipeline.canonical_stage("qc_review") == "post_qc_review"
    assert pipeline.canonical_stage("s1_rna_qc") == "s1_rna_qc"   # passthrough
    assert pipeline.display_stage("post_qc_review") == "qc_review"
    assert pipeline.display_stage("resolution_review") == "resolution_review"  # passthrough


def test_stages_follow_workflow_order():
    assert pipeline.PLANNING_STAGES == ("p1_context", "s0_ingest", "plan_review")
    assert pipeline.STAGES.index("post_qc_review") < pipeline.STAGES.index("qc_handoff")
    assert pipeline.STAGES.index("qc_handoff") < pipeline.STAGES.index("s4_rna_norm")


def test_branch_has_atac_reexport():
    assert pipeline.branch_has_atac("paired")
    assert pipeline.branch_has_atac("atac_only")
    assert not pipeline.branch_has_atac("rna_only")
