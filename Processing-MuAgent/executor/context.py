"""Biological Context Report parser + structured field extraction.

Non-LLM MVP: parses the markdown template, applies deterministic inference rules
(e.g. 'snRNA' -> sample_type=nuclei), labels each field with source + confidence + status.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from . import doi_fetch

TEMPLATE = """Biological Context Report
- Organism:
- Tissue / sample:
- Assay:
- DOI(s) of related or original paper:
- Anything else relevant (optional):
"""


def write_template(path: Path | str) -> None:
    p = Path(path)
    if p.exists():
        return
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(TEMPLATE)


def is_unfilled_template(path: Path | str) -> bool:
    """True if the Biological Context Report exists but has no user-provided values.

    A line like '- Organism:' (nothing after the colon) counts as unfilled. If every
    expected field line ends with a colon (no value), the template is considered
    unfilled and the subagent must stop and ask the user.
    """
    p = Path(path)
    if not p.exists():
        return True
    text = p.read_text()
    parsed = parse_report(text)
    # Consider it "unfilled" only if ALL four non-optional fields are empty.
    required = ["organism", "tissue", "assay"]
    # dois is optional — do not block on it.
    any_filled = any(parsed.get(k) for k in required)
    return not any_filled


class ContextNotProvided(RuntimeError):
    """Raised when the Biological Context Report is absent or still the blank template."""


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

_FIELD_MAP = {
    # Use [ \t]* (horizontal whitespace only) around the colon — `\s*` would
    # otherwise consume newlines and leak the next line's content into the
    # value when the field is blank.
    "organism": re.compile(r"^-[ \t]*Organism[ \t]*:[ \t]*(.*)$", re.IGNORECASE | re.MULTILINE),
    "tissue": re.compile(r"^-[ \t]*Tissue.*?:[ \t]*(.*)$", re.IGNORECASE | re.MULTILINE),
    "assay": re.compile(r"^-[ \t]*Assay[ \t]*:[ \t]*(.*)$", re.IGNORECASE | re.MULTILINE),
    "dois_raw": re.compile(r"^-[ \t]*DOI.*?:[ \t]*(.*)$", re.IGNORECASE | re.MULTILINE),
    "other": re.compile(r"^-[ \t]*Anything else.*?:[ \t]*(.*)$", re.IGNORECASE | re.MULTILINE),
}


def parse_report(text: str) -> dict[str, str | list[str]]:
    out: dict[str, str | list[str]] = {}
    for k, pat in _FIELD_MAP.items():
        m = pat.search(text)
        out[k] = m.group(1).strip() if m else ""
    # DOIs: split on whitespace / commas
    dois_raw = out.get("dois_raw", "") or ""
    dois = [d.strip() for d in re.split(r"[,\s]+", str(dois_raw)) if d.strip()]
    # Keep only tokens that look like DOIs (contain '10.' prefix)
    out["dois"] = [d for d in dois if "10." in d]
    return out


# ---------------------------------------------------------------------------
# Inference rules
# ---------------------------------------------------------------------------

def _infer_sample_type(assay: str) -> tuple[str, str, str]:
    """Return (value, confidence, rationale) for sample_type given assay string."""
    a = assay.lower()
    if "snrna" in a or "snatac" in a or "nuclei" in a or "single-nucleus" in a or "single nucleus" in a:
        return "nuclei", "high", "Assay string names single-nucleus protocol."
    if "scrna" in a or "scatac" in a or "single-cell" in a or "single cell" in a:
        return "cells", "high", "Assay string names single-cell protocol."
    return "unknown", "low", "Assay string does not indicate cells vs nuclei."


def _infer_modality_type(assay: str) -> tuple[str, str, str]:
    a = assay.lower()
    has_rna = any(t in a for t in ["rna", "gex", "gene expression"])
    has_atac = "atac" in a
    if has_rna and has_atac:
        return "paired_multiome_candidate", "medium", "Assay describes both RNA and ATAC; confirm pairing after ingest."
    if has_rna and not has_atac:
        return "rna_only", "high", "Assay names RNA only."
    if has_atac and not has_rna:
        return "atac_only", "high", "Assay names ATAC only."
    return "unknown", "low", "Assay does not clearly indicate modalities."


# ---------------------------------------------------------------------------
# Field-level extraction entry
# ---------------------------------------------------------------------------

def _entry(value: Any, source: str, confidence: str, status: str, rationale: str,
           source_refs: list[str] | None = None) -> dict[str, Any]:
    return {
        "value": value,
        "source": source,
        "confidence": confidence,
        "status": status,
        "rationale": rationale,
        "source_refs": source_refs or [],
    }


def extract_context(
    report_path: Path | str,
    run_dir: Path | str,
    *,
    file_input_signals: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Produce context_extraction.json payload."""
    report_path = Path(report_path)
    report_text = report_path.read_text() if report_path.exists() else ""
    parsed = parse_report(report_text)
    file_input_signals = file_input_signals or {}

    fields: dict[str, dict[str, Any]] = {}

    # Organism
    if parsed.get("organism"):
        fields["organism"] = _entry(
            parsed["organism"], "report", "high", "explicit",
            "User-declared in biological context.", [str(report_path)]
        )
    else:
        fields["organism"] = _entry("unknown", "inferred", "low", "missing",
                                    "Not declared; no reliable inference source.")

    # Tissue
    if parsed.get("tissue"):
        fields["tissue"] = _entry(parsed["tissue"], "report", "high", "explicit",
                                  "User-declared.", [str(report_path)])
    else:
        fields["tissue"] = _entry("unknown", "inferred", "low", "missing",
                                  "Not declared.")

    # Assay
    if parsed.get("assay"):
        fields["assay_type"] = _entry(parsed["assay"], "report", "high", "explicit",
                                      "User-declared.", [str(report_path)])
    else:
        fields["assay_type"] = _entry("unknown", "inferred", "low", "missing",
                                      "Not declared.")

    # Derived: sample_type (cells vs nuclei)
    assay_txt = str(parsed.get("assay", "") or "")
    st_value, st_conf, st_rationale = _infer_sample_type(assay_txt)
    fields["sample_type"] = _entry(
        st_value,
        "inferred" if st_value != "unknown" else "inferred",
        st_conf,
        "inferred" if st_value != "unknown" else "missing",
        st_rationale,
        [str(report_path)],
    )

    # Genome reference: user-declared in run.yaml only (no organism-based default).
    file_genome = file_input_signals.get("genome_assembly")
    if file_genome:
        fields["genome_build"] = _entry(
            file_genome,
            "run_yaml",
            "high",
            "explicit",
            "User-declared.",
        )
        gb_value = file_genome
    else:
        fields["genome_build"] = _entry(
            "unknown",
            "run_yaml",
            "low",
            "missing",
            "Declare the reference genome as part of the biological context.",
        )
        gb_value = "unknown"

    # Derived: modality_type
    mt_value, mt_conf, mt_rationale = _infer_modality_type(assay_txt)
    # Upgrade to file_input if signals say so
    if file_input_signals.get("both_modalities_in_single_h5"):
        mt_value = "paired_multiome"
        fields["modality_type"] = _entry(mt_value, "file_input", "high", "explicit",
                                          "Both Gene Expression and Peaks found in single Cell Ranger ARC .h5.",
                                          ["artifacts/s0_ingest/validation_report.json"])
    else:
        fields["modality_type"] = _entry(mt_value, "inferred", mt_conf,
                                         "inferred" if mt_value != "unknown" else "missing",
                                         mt_rationale, [str(report_path)])

    # DOI handling
    dois = parsed.get("dois") or []
    doi_entries: list[dict[str, Any]] = []
    for doi in dois:
        fetched = doi_fetch.fetch_and_cache(run_dir, doi)
        doi_entries.append({
            "doi": doi,
            "slug": doi_fetch.doi_slug(doi),
            "status": fetched.get("status"),
            "title": fetched.get("title"),
            "year": fetched.get("year"),
        })
    fields["dois"] = _entry(
        doi_entries,
        "report" if doi_entries else "inferred",
        "high" if doi_entries else "low",
        "explicit" if doi_entries else "missing",
        f"Fetched {len(doi_entries)} DOI(s) via Crossref." if doi_entries else "No DOIs supplied.",
        [str(report_path)] if doi_entries else [],
    )

    conflicts: list[dict[str, Any]] = []
    file_chroms = set(file_input_signals.get("fragment_chromosomes", []) or [])
    if file_chroms:
        has_chr = any(c.startswith("chr") for c in file_chroms)
        if gb_value == "mm10" and not has_chr and all(not c.startswith("chr") for c in file_chroms):
            # Ensembl-style naming — not a conflict per se, but note it
            pass

    return {
        "fields": fields,
        "conflicts": conflicts,
        "report_text": report_text,
        "parsed": parsed,
    }


def write_context_extraction(run_dir: Path | str, payload: dict[str, Any]) -> Path:
    from .run_paths import RunPaths
    out = RunPaths(Path(run_dir)).artifact("p1_context", "context_extraction.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        json.dump(payload, f, indent=2, default=str)
    return out


def render_summary(payload: dict[str, Any]) -> str:
    lines: list[str] = ["# Context Extraction Summary", ""]
    lines.append("## Fields")
    for k, v in payload["fields"].items():
        if k == "dois":
            continue
        val_str = v["value"] if not isinstance(v["value"], list) else ", ".join(map(str, v["value"]))
        lines.append(f"- **{k}**: `{val_str}` — {v['source']} / {v['confidence']} / {v['status']}")
        lines.append(f"  - {v['rationale']}")
    lines.append("")
    lines.append("## DOIs")
    dois = payload["fields"].get("dois", {}).get("value") or []
    if dois:
        for d in dois:
            lines.append(f"- {d['doi']} — {d.get('status')} — {d.get('title', '')}")
    else:
        lines.append("- (none)")
    lines.append("")
    lines.append("## Conflicts")
    if payload["conflicts"]:
        for c in payload["conflicts"]:
            vals = ", ".join(f"{v['source']}={v['value']}" for v in c["values"])
            lines.append(f"- **{c['field']}** [{c['severity']}]: {vals}")
    else:
        lines.append("- (none)")
    return "\n".join(lines) + "\n"
