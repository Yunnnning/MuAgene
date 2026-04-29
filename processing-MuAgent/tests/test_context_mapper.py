"""context_mapper round-trips through context.parse_report + write_report path."""
from __future__ import annotations

import tempfile
from pathlib import Path

from executor import context, context_mapper


def test_build_roundtrips_through_parse_report() -> None:
    md = context_mapper.build_report_from_chat(
        organism="mouse", tissue="testis",
        assay="single-nucleus multiome (snRNA + snATAC)",
        dois=["10.1016/j.stemcr.2025.102449"],
        notes="GSE268104",
    )
    parsed = context.parse_report(md)
    assert parsed["organism"] == "mouse"
    assert parsed["tissue"] == "testis"
    assert "snRNA" in parsed["assay"]
    assert parsed["dois"] == ["10.1016/j.stemcr.2025.102449"]


def test_filled_report_is_not_unfilled() -> None:
    md = context_mapper.build_report_from_chat(
        organism="mouse", tissue="testis", assay="snRNA-seq",
    )
    with tempfile.NamedTemporaryFile(suffix=".md", mode="w", delete=False) as f:
        f.write(md)
        p = f.name
    try:
        assert context.is_unfilled_template(p) is False
    finally:
        Path(p).unlink(missing_ok=True)


def test_canonical_template_is_unfilled() -> None:
    """Regression guard for the `\\s*`→`[ \\t]*` regex fix in `parse_report`."""
    with tempfile.NamedTemporaryFile(suffix=".md", mode="w", delete=False) as f:
        f.write(context.TEMPLATE)
        p = f.name
    try:
        assert context.is_unfilled_template(p) is True
    finally:
        Path(p).unlink(missing_ok=True)


def test_empty_build_is_unfilled() -> None:
    empty = context_mapper.build_report_from_chat()
    with tempfile.NamedTemporaryFile(suffix=".md", mode="w", delete=False) as f:
        f.write(empty)
        p = f.name
    try:
        assert context.is_unfilled_template(p) is True
    finally:
        Path(p).unlink(missing_ok=True)


def test_append_dois_dedup() -> None:
    md = context_mapper.build_report_from_chat(
        organism="mouse", tissue="testis", assay="snRNA",
        dois=["10.1016/j.stemcr.2025.102449"],
    )
    md2 = context_mapper.append_dois(md, ["10.9999/foo", "10.1016/j.stemcr.2025.102449"])
    parsed = context.parse_report(md2)
    assert set(parsed["dois"]) == {"10.1016/j.stemcr.2025.102449", "10.9999/foo"}


def test_write_report_lands_in_pre_run_config(tmp_path: Path) -> None:
    md = context_mapper.build_report_from_chat(organism="mouse", tissue="testis", assay="scRNA")
    out = context_mapper.write_report(tmp_path, md)
    assert out == tmp_path / "deliverables" / "pre_run" / "config" / "biological_context.md"
    assert out.read_text().startswith("Biological Context Report")
