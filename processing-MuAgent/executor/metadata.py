"""Metadata ingestion: provided-file validation, recovery, minimal reconstruction."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import pandas as pd


def load_user_metadata(path: Path | str) -> pd.DataFrame:
    p = Path(path)
    sep = "\t" if p.suffix in {".tsv", ".txt"} else ","
    return pd.read_csv(p, sep=sep)


def identify_join_key(df: pd.DataFrame, rna_bc: set[str], atac_bc: set[str]) -> tuple[str, float]:
    """Pick the column whose values overlap barcode sets most."""
    best_col, best_cov = None, 0.0
    n_total = max(len(rna_bc | atac_bc), 1)
    for col in df.columns:
        vals = set(df[col].astype(str))
        cov = len(vals & (rna_bc | atac_bc)) / n_total
        if cov > best_cov:
            best_cov, best_col = cov, col
    return (best_col or df.columns[0]), best_cov


def reconstruct_minimal(
    rna_bc: set[str],
    atac_bc: set[str],
    out_path: Path | str,
    *,
    inferred_sample_id_fn=None,
) -> pd.DataFrame:
    """Build metadata_minimal.tsv with cell_id, modality, inferred_sample_id, batch, condition."""
    union = sorted(rna_bc | atac_bc)
    rows: list[dict[str, Any]] = []
    for bc in union:
        in_rna = bc in rna_bc
        in_atac = bc in atac_bc
        modality = "both" if (in_rna and in_atac) else ("rna" if in_rna else "atac")
        sid = inferred_sample_id_fn(bc) if inferred_sample_id_fn else ""
        rows.append({
            "cell_id": bc,
            "modality": modality,
            "inferred_sample_id": sid,
            "batch": "unknown",
            "condition": "unknown",
        })
    df = pd.DataFrame(rows)
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    return df


def search_adjacent_metadata(input_dirs: list[Path]) -> list[dict[str, Any]]:
    """Look for sibling metadata-ish files; return list of candidates with path + reason."""
    hits: list[dict[str, Any]] = []
    patterns = ["*metadata*", "*samples*", "GSE*_series_matrix*", "*.csv", "*.tsv"]
    for d in input_dirs:
        if not d.exists() or not d.is_dir():
            continue
        for pat in patterns:
            for f in d.glob(pat):
                if f.is_file() and f.stat().st_size > 0:
                    hits.append({"path": str(f), "reason": f"glob:{pat}"})
    return hits


def unrecoverable_categories(source: str) -> list[str]:
    """Categories that require user-supplied metadata; always listed to surface limitations."""
    if source == "provided":
        return []
    return ["batch", "condition", "donor_grouping"]
