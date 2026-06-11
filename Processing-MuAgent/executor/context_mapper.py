"""Build a Biological Context Report from already-extracted fields.

The chat-agent runtime owns natural-language understanding: it asks the user,
interprets free text, and hands the executor a structured tuple of fields. This
module is a deterministic markdown emitter — no LLM calls, no heuristic parsing.
The output format matches `executor.context.TEMPLATE` so it round-trips cleanly
through `context.parse_report` and satisfies `context.is_unfilled_template`.

Typical agent use:

    from executor import context_mapper
    md = context_mapper.build_report_from_chat(
        organism="mouse", tissue="testis",
        assay="single-nucleus multiome (snRNA + snATAC)",
        dois=["10.1016/j.stemcr.2025.102449"],
        notes="GSE268104, adult mouse C57BL/6")
    context_mapper.write_report(run_dir, md)
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable


def build_report_from_chat(
    *,
    organism: str = "",
    tissue: str = "",
    assay: str = "",
    dois: Iterable[str] | None = None,
    notes: str = "",
) -> str:
    """Return a Biological Context Report markdown string.

    All arguments are keyword-only — call sites must name them so it is obvious
    at the call point which field each value populates. Empty strings emit an
    empty field (same shape as the blank template). DOI strings are joined
    with commas; non-DOI tokens are allowed but `context.parse_report` filters
    them out when it reads the file back.
    """
    doi_line = ", ".join(d.strip() for d in (dois or []) if d.strip())
    lines = [
        "Biological Context Report",
        f"- Organism: {organism.strip()}",
        f"- Tissue / sample: {tissue.strip()}",
        f"- Assay: {assay.strip()}",
        f"- DOI(s) of related or original paper: {doi_line}",
        f"- Anything else relevant (optional): {notes.strip()}",
    ]
    return "\n".join(lines) + "\n"


def append_dois(report_md: str, dois: Iterable[str]) -> str:
    """Merge additional DOI tokens into an existing report markdown string.

    Preserves existing DOIs (order + formatting); appends new ones separated by
    commas. If the existing DOI line is empty, writes the new list verbatim.
    Unknown lines pass through unchanged.
    """
    new_tokens = [d.strip() for d in dois if d.strip()]
    if not new_tokens:
        return report_md

    out_lines: list[str] = []
    touched = False
    for line in report_md.splitlines():
        stripped = line.strip()
        lower = stripped.lower()
        if not touched and lower.startswith("- doi"):
            # Line shape: "- DOI(s) of related or original paper: <existing>"
            head, _, existing = line.partition(":")
            existing_tokens = [d.strip() for d in existing.split(",") if d.strip()]
            merged = existing_tokens + [d for d in new_tokens if d not in existing_tokens]
            out_lines.append(f"{head}: {', '.join(merged)}")
            touched = True
        else:
            out_lines.append(line)
    return "\n".join(out_lines) + ("\n" if report_md.endswith("\n") else "")


def write_report(run_dir: Path | str, content: str) -> Path:
    """Write `content` to the canonical Biological Context Report path.

    Resolves to `RunPaths(run_dir).biological_context_md` — currently
    `<run_dir>/deliverables/plan/config/biological_context.md`. Creates
    parent directories as needed. Returns the written path.
    """
    from .run_paths import RunPaths
    out = RunPaths(Path(run_dir)).biological_context_md
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(content)
    return out
