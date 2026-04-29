"""Branch-aware wording for the doublet-policy review item.

Contract:
- `paired` / `separate` → certainty "needs confirmation"; wording says
  reconciliation.
- `rna_only` / `atac_only` → certainty "certain"; wording says auto-applied,
  single-detector; no confirmation implied.
"""
from __future__ import annotations

import json
from pathlib import Path

from executor.plan_review import build_summary


def _fixture_run(tmp_path: Path, *, branch: str) -> Path:
    """Seed the minimal artifacts build_summary reads: p1 context, s0 ingest, p2 plan."""
    run_dir = tmp_path
    art = run_dir / "internal" / "artifacts"
    (art / "p1_context").mkdir(parents=True, exist_ok=True)
    (art / "s0_ingest").mkdir(parents=True, exist_ok=True)
    (art / "p2_plan").mkdir(parents=True, exist_ok=True)

    (art / "p1_context" / "context_extraction.json").write_text(json.dumps(
        {"fields": {
            "modality_type": {"value": branch, "rationale": "test", "status": "explicit"},
            "organism": {"value": "mouse", "rationale": "test", "status": "explicit"},
            "genome_build": {"value": "mm10", "rationale": "test", "confidence": "high"},
            "sample_type": {"value": "nuclei"},
        }}
    ))
    (art / "s0_ingest" / "validation_report.json").write_text(json.dumps(
        {"pairing": {"status": branch, "overlap": 0.0, "method": "test",
                     "confidence": "high", "n_shared": 0}}
    ))
    (art / "p2_plan" / "preprocessing_plan.json").write_text(json.dumps({
        "workflow_branch": branch,
        "stages": {
            "s1_rna_qc": {"parameters": {"pct_mt_ceiling": {"value": 10}}},
            "s2_atac_qc": {"parameters": {"tss_enrichment_min": {"value": 2}}},
            "s3_doublets": {"parameters": {
                "removal_policy_recommendation": {"value": "union"},
                "study_goal": {"value": "clustering_inference"},
            }},
            "s7_clustering": {"parameters": {
                "leiden_resolution_grid": {"value": [0.6, 0.8, 1.0]},
                "rna_tilt": {"value": "higher"},
                "atac_tilt": {"value": "lower"},
            }},
        },
    }))
    return run_dir


def _doublet_item(items: list[dict]) -> dict:
    [item] = [i for i in items if i["label"] == "Doublet removal policy"]
    return item


def test_paired_needs_confirmation(tmp_path: Path) -> None:
    items = build_summary(_fixture_run(tmp_path, branch="paired"))
    item = _doublet_item(items)
    assert item["certainty"] == "needs confirmation"
    assert "reconcil" in item["value"].lower() or "reconcil" in item["reason"].lower()


def test_separate_needs_confirmation(tmp_path: Path) -> None:
    items = build_summary(_fixture_run(tmp_path, branch="separate"))
    item = _doublet_item(items)
    assert item["certainty"] == "needs confirmation"


def test_rna_only_auto_applied(tmp_path: Path) -> None:
    items = build_summary(_fixture_run(tmp_path, branch="rna_only"))
    item = _doublet_item(items)
    assert item["certainty"] == "certain"
    assert "auto-applied" in item["value"].lower()
    assert "scrublet" in item["value"].lower()
    # must not suggest user confirmation is required
    assert "needs confirmation" not in item["certainty"]


def test_atac_only_auto_applied(tmp_path: Path) -> None:
    items = build_summary(_fixture_run(tmp_path, branch="atac_only"))
    item = _doublet_item(items)
    assert item["certainty"] == "certain"
    assert "auto-applied" in item["value"].lower()
    assert "atac" in item["value"].lower()


def test_no_keep_all_wording_in_any_branch(tmp_path: Path) -> None:
    """Regression guard: `keep-all` is not implemented end-to-end and must not leak into user-facing wording."""
    for branch in ("paired", "separate", "rna_only", "atac_only"):
        items = build_summary(_fixture_run(tmp_path, branch=branch))
        item = _doublet_item(items)
        combined = f"{item['value']} {item['reason']}".lower()
        assert "keep-all" not in combined and "keep_all" not in combined, (
            f"keep-all leaked into doublet-policy wording for branch={branch}: {item}"
        )
