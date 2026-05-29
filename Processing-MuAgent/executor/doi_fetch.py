"""DOI metadata fetch via Crossref (polite API, no auth required).

Classifies each DOI as fetched | abstract_only | fetch_failed | paywalled.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import requests


CROSSREF_URL = "https://api.crossref.org/works/{doi}"
HEADERS = {"User-Agent": "Processing-MuAgent/0.1 (mailto:yuanyunning@gmail.com)"}


def doi_slug(doi: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", doi.strip())


def fetch_doi(doi: str, timeout: float = 15.0) -> dict[str, Any]:
    try:
        r = requests.get(CROSSREF_URL.format(doi=doi), headers=HEADERS, timeout=timeout)
    except Exception as e:
        return {"doi": doi, "status": "fetch_failed", "error": str(e)}
    if r.status_code == 404:
        return {"doi": doi, "status": "fetch_failed", "error": "DOI not in Crossref"}
    if r.status_code != 200:
        return {"doi": doi, "status": "fetch_failed", "error": f"HTTP {r.status_code}"}
    try:
        data = r.json().get("message", {})
    except Exception:
        return {"doi": doi, "status": "fetch_failed", "error": "non-JSON response"}
    abstract = data.get("abstract")  # Crossref returns JATS-tagged string when present
    title = (data.get("title") or [""])[0]
    authors = [f"{a.get('family', '')}, {a.get('given', '')}" for a in data.get("author", [])]
    journal = (data.get("container-title") or [""])[0]
    year = None
    for key in ("published-print", "published-online", "published"):
        if key in data:
            parts = data[key].get("date-parts", [[None]])[0]
            year = parts[0] if parts else None
            if year:
                break
    return {
        "doi": doi,
        "status": "fetched" if abstract else "abstract_only" if title else "fetch_failed",
        "title": title,
        "authors": authors,
        "journal": journal,
        "year": year,
        "abstract": abstract,
        "raw": data,
    }


def cache_path(run_dir: Path | str, doi: str) -> Path:
    from .run_paths import RunPaths
    return RunPaths(Path(run_dir)).artifact("p1_context", f"doi_{doi_slug(doi)}.json")


def fetch_and_cache(run_dir: Path | str, doi: str) -> dict[str, Any]:
    out = cache_path(run_dir, doi)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        with out.open() as f:
            return json.load(f)
    result = fetch_doi(doi)
    with out.open("w") as f:
        json.dump(result, f, indent=2)
    return result
