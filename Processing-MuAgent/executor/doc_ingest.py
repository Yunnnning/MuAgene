"""Ingest a user-supplied Biological Context Report document.

Supported formats: .docx, .pdf, .md, .txt.

Returns plain text; downstream parsing (context.parse_report) is format-agnostic
as long as the document contains the standard template lines (`- Organism:`, etc.)
or can be matched heuristically.
"""
from __future__ import annotations

import re
from pathlib import Path


def _read_docx(path: Path) -> str:
    import docx  # python-docx
    d = docx.Document(str(path))
    lines: list[str] = []
    for para in d.paragraphs:
        t = para.text.strip()
        if t:
            lines.append(t)
    # Tables (occasionally used for the template)
    for table in getattr(d, "tables", []):
        for row in table.rows:
            cells = [c.text.strip() for c in row.cells if c.text.strip()]
            if cells:
                lines.append(" | ".join(cells))
    return "\n".join(lines)


def _read_pdf(path: Path) -> str:
    from pypdf import PdfReader
    r = PdfReader(str(path))
    out: list[str] = []
    for page in r.pages:
        t = page.extract_text()
        if t:
            out.append(t)
    return "\n".join(out)


def _read_text(path: Path) -> str:
    return Path(path).read_text(encoding="utf-8", errors="replace")


def read_document(path: Path | str) -> str:
    """Return plain text from a supported document. Raises on unsupported extension."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    suffix = p.suffix.lower()
    if suffix == ".docx":
        return _read_docx(p)
    if suffix == ".pdf":
        return _read_pdf(p)
    if suffix in {".md", ".txt"}:
        return _read_text(p)
    raise ValueError(f"Unsupported document format: {suffix}. Accepted: .docx, .pdf, .md, .txt")


# ---------------------------------------------------------------------------
# Heuristic field extraction — the document may or may not match the strict
# template. We run the strict regex first, then a looser label-matching pass
# for paragraphs like "Organism: mouse." (period-terminated, no leading dash).
# ---------------------------------------------------------------------------

_LABEL_PATTERNS = {
    "organism":  re.compile(r"(?:^|\n)\s*(?:-\s*)?organism\s*:\s*(.+?)(?:\n|\.$|$)", re.IGNORECASE),
    "tissue":    re.compile(r"(?:^|\n)\s*(?:-\s*)?tissue[^:]*:\s*(.+?)(?:\n|\.$|$)", re.IGNORECASE),
    "assay":     re.compile(r"(?:^|\n)\s*(?:-\s*)?assay\s*:\s*(.+?)(?:\n|\.$|$)", re.IGNORECASE),
    "doi_line":  re.compile(r"(?:^|\n)\s*(?:-\s*)?DOI(?:\(s\))?[^:]*:\s*(.+?)(?:\n|$)", re.IGNORECASE),
    "notes":     re.compile(r"(?:^|\n)\s*(?:-\s*)?(?:additional\s*notes|anything\s*else[^:]*|notes)\s*:\s*(.+?)(?:\n|$)", re.IGNORECASE),
}

_DOI_REGEX = re.compile(r"(10\.\d{4,9}/[^\s,;]+)", re.IGNORECASE)


def extract_fields(text: str) -> dict[str, object]:
    """Return a best-effort dict of fields extracted from plain text."""
    out: dict[str, object] = {"organism": "", "tissue": "", "assay": "",
                               "dois": [], "notes": "", "dois_raw": ""}
    for key, pat in _LABEL_PATTERNS.items():
        m = pat.search(text)
        if not m:
            continue
        val = m.group(1).strip().rstrip(".,;")
        if key == "doi_line":
            out["dois_raw"] = val
            dois = _DOI_REGEX.findall(val)
            if not dois:
                # The value may itself be a URL that embeds a DOI
                dois = _DOI_REGEX.findall(val.replace("https://doi.org/", ""))
            out["dois"] = dois
        elif key == "notes":
            out["notes"] = val
        else:
            out[key] = val
    # Fallback DOI scan on whole text if the labelled line missed
    if not out["dois"]:
        out["dois"] = list(set(_DOI_REGEX.findall(text)))
    return out


def canonicalise_to_template(fields: dict[str, object]) -> str:
    """Render extracted fields as the canonical biological_context.md template.

    Used to write a reconstructed biological_context.md beside the original doc
    so downstream stages have a single source of truth.
    """
    dois_str = ", ".join(fields.get("dois") or []) or str(fields.get("dois_raw") or "")
    notes = fields.get("notes") or ""
    return (
        "Biological Context Report\n"
        f"- Organism: {fields.get('organism', '')}\n"
        f"- Tissue / sample: {fields.get('tissue', '')}\n"
        f"- Assay: {fields.get('assay', '')}\n"
        f"- DOI(s) of related or original paper: {dois_str}\n"
        f"- Anything else relevant (optional): {notes}\n"
    )
