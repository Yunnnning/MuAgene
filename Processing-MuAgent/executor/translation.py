"""Barcode translation between ATAC and RNA whitelist spaces.

Paired multiome runs produced by separate Cell Ranger (GEX) + Cell Ranger
ATAC pipelines end up with cell barcodes in two different 10x whitelists
(GEX whitelist vs ATAC whitelist), which never overlap algorithmically. When
the user supplies a translation table — a 2-column TSV mapping
`rna_barcode <-> atac_barcode` per cell — S0's diagnostics ladder uses it to
rewrite the ATAC fragment barcodes into RNA-space, so that downstream
intersection at S3 actually finds shared cells.

The translation table itself is the *only* artifact this module writes that
later stages consult: it lives at
`<run_dir>/internal/artifacts/s0_ingest/barcode_translation.parquet` and is
read by S2 to produce a translated copy of `atac_fragments.tsv.gz` upstream
of the SnapATAC2 import call.

S2's QC code path is byte-identical regardless of translation — only the
file it reads differs.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd


def load_translation_tsv(path: Path | str) -> dict[str, str]:
    """Read a 2-column TSV and return the ATAC->RNA mapping.

    Accepted column-name conventions (case-insensitive):
      - `rna_barcode`, `atac_barcode`
      - `gex_barcode`, `atac_barcode`
      - first two columns regardless of name (last resort).

    Rows with empty / missing values on either side are dropped. Duplicate
    ATAC keys keep the first occurrence (with a warning logged by the
    caller when count > 0).
    """
    p = Path(path)
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    if df.shape[1] < 2:
        raise ValueError(
            f"translation: {p} must have at least 2 columns (rna_barcode, atac_barcode); "
            f"got {df.shape[1]}."
        )
    cols_lower = {c.lower(): c for c in df.columns}
    rna_col = cols_lower.get("rna_barcode") or cols_lower.get("gex_barcode")
    atac_col = cols_lower.get("atac_barcode")
    if rna_col is None or atac_col is None:
        # Fall back to first two columns.
        rna_col, atac_col = df.columns[0], df.columns[1]
    pairs = df[[atac_col, rna_col]].astype(str)
    pairs = pairs[(pairs[atac_col] != "") & (pairs[rna_col] != "")]
    out: dict[str, str] = {}
    for atac_bc, rna_bc in zip(pairs[atac_col], pairs[rna_col]):
        out.setdefault(atac_bc, rna_bc)
    return out


def write_translation_parquet(table: dict[str, str], path: Path | str) -> Path:
    """Persist the ATAC->RNA mapping as parquet for fast re-read by S2."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        {"atac_barcode": list(table.keys()), "rna_barcode": list(table.values())}
    )
    df.to_parquet(p, index=False)
    return p


def load_translation_parquet(path: Path | str) -> dict[str, str]:
    df = pd.read_parquet(path)
    return dict(zip(df["atac_barcode"].astype(str), df["rna_barcode"].astype(str)))


def translate_atac_barcode_set(
    atac_bc: set[str], table: dict[str, str]
) -> tuple[set[str], int]:
    """Apply ATAC->RNA mapping to a set of ATAC barcodes.

    Returns (translated_set_in_rna_space, n_unmapped). Unmapped ATAC barcodes
    are silently dropped from the returned set; the count is surfaced so S0
    can record coverage in `validation_report.json`.
    """
    translated: set[str] = set()
    n_unmapped = 0
    for bc in atac_bc:
        rna = table.get(bc)
        if rna is None:
            n_unmapped += 1
        else:
            translated.add(rna)
    return translated, n_unmapped


def translate_fragments_file(
    src_path: Path | str,
    dst_path: Path | str,
    table: dict[str, str],
) -> dict[str, int]:
    """Stream-rewrite ATAC fragments, mapping column 4 (barcode) into RNA-space.

    Fragments whose ATAC barcode is not in `table` are dropped. The output
    file is gzip-compressed at default level, matching the input convention.
    Returns `{n_in, n_out, n_dropped}`. No tabix index is produced — SnapATAC2's
    `pp.import_fragments` reads sequentially and does not require one.
    """
    src_path = Path(src_path)
    dst_path = Path(dst_path)
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    n_in = 0
    n_out = 0
    with gzip.open(src_path, "rt") as fin, gzip.open(dst_path, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            n_in += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            mapped = table.get(parts[3])
            if mapped is None:
                continue
            parts[3] = mapped
            fout.write("\t".join(parts) + "\n")
            n_out += 1
    return {"n_in": n_in, "n_out": n_out, "n_dropped": n_in - n_out}
