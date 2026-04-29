"""ATAC-QC consistency guards.

Three contracts:
  1. Plan no longer exposes `nucleosome_signal_max` (metric is not reliably
     implemented in the SnapATAC2 path used here).
  2. plan_review QC-strategy wording does not claim a nucleosome-signal filter.
  3. s2_atac_qc rule proposal text does not claim a tile-matrix step (tile
     construction lives in S5).
"""
from __future__ import annotations

import json
from pathlib import Path

from executor.plan_assembler import assemble_plan
from executor.plan_review import build_summary


RULES = Path(__file__).parent.parent / "workflow" / "rules"


def test_plan_has_no_nucleosome_signal_max() -> None:
    for branch in ("paired", "separate", "atac_only"):
        plan = assemble_plan("/tmp/pma_noop", workflow_branch=branch)
        params = plan["stages"]["s2_atac_qc"]["parameters"]
        assert "nucleosome_signal_max" not in params, (
            f"nucleosome_signal_max leaked into the {branch} plan; remove it — "
            "the metric is not reliably computed by the current SnapATAC2 path."
        )


def test_plan_review_qc_wording_no_nucleosome(tmp_path: Path) -> None:
    art = tmp_path / "internal" / "artifacts"
    (art / "p1_context").mkdir(parents=True)
    (art / "s0_ingest").mkdir(parents=True)
    (art / "p2_plan").mkdir(parents=True)
    (art / "p1_context" / "context_extraction.json").write_text(json.dumps({
        "fields": {
            "modality_type": {"value": "paired", "rationale": "", "status": "explicit"},
            "organism": {"value": "mouse", "status": "explicit"},
            "genome_build": {"value": "mm10", "confidence": "high"},
            "sample_type": {"value": "nuclei"},
        }
    }))
    (art / "s0_ingest" / "validation_report.json").write_text(json.dumps({
        "pairing": {"status": "paired", "overlap": 0.99, "method": "exact",
                     "confidence": "high", "n_shared": 100}
    }))
    plan = assemble_plan(tmp_path, workflow_branch="paired", sample_type="nuclei")
    (art / "p2_plan" / "preprocessing_plan.json").write_text(json.dumps(plan))

    items = build_summary(tmp_path)
    qc = next(i for i in items if i["label"] == "QC strategy")
    assert "nucleosome" not in qc["value"].lower()
    assert "nucleosome" not in qc["reason"].lower()


def test_s2_rule_proposal_wording_no_tile_matrix() -> None:
    text = (RULES / "s2_atac_qc.smk").read_text()
    # The propose rule's action string must not claim tile-matrix creation —
    # tiles are built in S5. "(no tile matrix here — S5 builds it)" is
    # explicitly allowed as disambiguation wording.
    lines = [line for line in text.splitlines() if '"action"' in line]
    assert lines, "expected an action string in s2_atac_qc.smk"
    for line in lines:
        lower = line.lower()
        if "tile" in lower:
            assert "no tile matrix here" in lower, (
                f"s2_atac_qc proposal action still claims tile-matrix work: {line!r}"
            )
