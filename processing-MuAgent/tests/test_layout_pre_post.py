"""RunPaths pre_run/post_run layout contract."""
from __future__ import annotations

from pathlib import Path

from executor.run_paths import RunPaths


def test_ensure_creates_pre_and_post_subtrees(tmp_path: Path) -> None:
    p = RunPaths(tmp_path)
    p.ensure()
    # Expected deliverable sub-directories
    assert (tmp_path / "deliverables" / "pre_run" / "config").is_dir()
    assert (tmp_path / "deliverables" / "pre_run" / "summary").is_dir()
    assert (tmp_path / "deliverables" / "post_run" / "summary").is_dir()
    assert (tmp_path / "deliverables" / "post_run" / "figures").is_dir()
    assert (tmp_path / "deliverables" / "post_run" / "processed").is_dir()
    assert (tmp_path / "deliverables" / "post_run" / "notebooks").is_dir()
    # Legacy top-level dirs must NOT be created
    for legacy in ("config", "summary", "figures", "processed", "notebooks"):
        assert not (tmp_path / "deliverables" / legacy).exists()


def test_canonical_paths_split_by_phase(tmp_path: Path) -> None:
    p = RunPaths(tmp_path)
    # pre_run
    assert p.run_yaml == tmp_path / "deliverables" / "pre_run" / "config" / "run.yaml"
    assert p.biological_context_md == tmp_path / "deliverables" / "pre_run" / "config" / "biological_context.md"
    assert p.context_summary_md == tmp_path / "deliverables" / "pre_run" / "summary" / "context_summary.md"
    assert p.plan_summary_md == tmp_path / "deliverables" / "pre_run" / "summary" / "plan_summary.md"
    assert p.plan_review_md == tmp_path / "deliverables" / "pre_run" / "summary" / "plan_review.md"
    # post_run
    assert p.resolution_summary_md == tmp_path / "deliverables" / "post_run" / "summary" / "resolution_summary.md"
    assert p.qc_summary_md == tmp_path / "deliverables" / "post_run" / "summary" / "qc_summary.md"
    assert p.run_manifest_json == tmp_path / "deliverables" / "post_run" / "summary" / "run_manifest.json"
    assert p.layout_json == tmp_path / "deliverables" / "post_run" / "summary" / "layout.json"
    assert p.processed_h5mu == tmp_path / "deliverables" / "post_run" / "processed" / "processed.h5mu"
    assert p.rna_processed_h5ad == tmp_path / "deliverables" / "post_run" / "processed" / "rna_processed.h5ad"
    assert p.atac_processed_h5ad == tmp_path / "deliverables" / "post_run" / "processed" / "atac_processed.h5ad"
