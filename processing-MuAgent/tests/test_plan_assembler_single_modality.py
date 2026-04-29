"""plan_assembler emits the correct per-branch stage set."""
from __future__ import annotations

import pytest

from executor.plan_assembler import _stages_for_branch, assemble_plan


def test_rna_only_stages() -> None:
    assert _stages_for_branch("rna_only") == {
        "s1_rna_qc", "s3_doublets", "s4_rna_norm",
        "s6_dimred", "s7_clustering", "s8_umap",
    }


def test_atac_only_stages() -> None:
    assert _stages_for_branch("atac_only") == {
        "s2_atac_qc", "s3_doublets", "s5_atac_lsi",
        "s6_dimred", "s7_clustering", "s8_umap",
    }


def test_paired_and_separate_stages_unchanged() -> None:
    full = {
        "s1_rna_qc", "s2_atac_qc", "s3_doublets", "s4_rna_norm",
        "s5_atac_lsi", "s6_dimred", "s7_clustering", "s8_umap",
    }
    assert _stages_for_branch("paired") == full
    assert _stages_for_branch("separate") == full


def test_unknown_branch_raises() -> None:
    with pytest.raises(ValueError, match="Unknown workflow_branch"):
        _stages_for_branch("wnn_integration")


def test_assemble_plan_filters_stages(tmp_path) -> None:
    plan = assemble_plan(tmp_path, workflow_branch="rna_only",
                         sample_type="nuclei", study_goal="clustering_inference")
    assert plan["workflow_branch"] == "rna_only"
    assert set(plan["stages"].keys()) == _stages_for_branch("rna_only")
    # Spot-check that a plan entry retains its provenance structure.
    s1 = plan["stages"]["s1_rna_qc"]["parameters"]
    assert "pct_mt_ceiling" in s1
    v = s1["pct_mt_ceiling"]
    assert set(v.keys()) >= {"value", "source", "rationale", "confidence"}


def test_assemble_plan_atac_only(tmp_path) -> None:
    plan = assemble_plan(tmp_path, workflow_branch="atac_only")
    assert plan["workflow_branch"] == "atac_only"
    assert set(plan["stages"].keys()) == _stages_for_branch("atac_only")
    # RNA-only stages must NOT appear
    assert "s1_rna_qc" not in plan["stages"]
    assert "s4_rna_norm" not in plan["stages"]
