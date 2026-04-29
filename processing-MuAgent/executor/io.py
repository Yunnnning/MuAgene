"""RNA/ATAC input format autodetection and loading.

Supported RNA formats: 10x HDF5 (.h5), 10x MEX dir, AnnData (.h5ad), custom.
ATAC: fragments.tsv.gz (+ .tbi required), optional peak matrix / h5ad.
"""
from __future__ import annotations

import gzip
import json
import subprocess
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import scanpy as sc


# ---------------------------------------------------------------------------
# RNA format autodetect
# ---------------------------------------------------------------------------

def detect_rna_format(path: Path | str) -> str:
    p = Path(path)
    if p.is_dir():
        names = {x.name for x in p.iterdir()}
        stems = {n.split(".")[0] for n in names}
        if ("matrix.mtx" in names or "matrix.mtx.gz" in names) and (
            "barcodes.tsv" in names or "barcodes.tsv.gz" in names
        ):
            return "10x_mex"
        raise ValueError(f"Directory {p} doesn't look like a 10x MEX bundle. Found: {names}")
    if p.suffix == ".h5ad":
        return "h5ad"
    if p.suffix == ".h5":
        # Peek inside: 10x Cell Ranger layout has /matrix group.
        try:
            with h5py.File(p, "r") as f:
                if "matrix" in f:
                    return "10x_h5"
        except Exception:
            pass
        raise ValueError(f"{p}: .h5 file lacks 10x Cell Ranger /matrix group")
    raise ValueError(f"Cannot autodetect RNA format for {p}")


def load_rna(path: Path | str, fmt: str | None = None) -> ad.AnnData:
    fmt = fmt or detect_rna_format(path)
    if fmt == "10x_h5":
        # Cell Ranger ARC h5 contains both Gene Expression and Peaks feature types.
        # sc.read_10x_h5 returns one AnnData with all features; we filter to Gene Expression.
        a = sc.read_10x_h5(str(path), gex_only=False)
        if "feature_types" in a.var.columns:
            a = a[:, a.var["feature_types"] == "Gene Expression"].copy()
        a.var_names_make_unique()
        return a
    if fmt == "10x_mex":
        a = sc.read_10x_mtx(str(path), var_names="gene_symbols", cache=False)
        a.var_names_make_unique()
        return a
    if fmt == "h5ad":
        return ad.read_h5ad(str(path))
    raise ValueError(f"Unknown RNA format {fmt!r}")


def detect_peaks_in_10x_h5(path: Path | str) -> bool:
    """Return True if the 10x h5 has 'Peaks' feature type (Cell Ranger ARC)."""
    try:
        with h5py.File(path, "r") as f:
            if "matrix" not in f:
                return False
            m = f["matrix"]
            if "features" not in m:
                return False
            if "feature_type" in m["features"]:
                ft = m["features"]["feature_type"][...]
                # bytes or str
                ft_set = {x.decode() if isinstance(x, bytes) else x for x in ft}
                return "Peaks" in ft_set
    except Exception:
        pass
    return False


def load_atac_from_10x_h5(path: Path | str) -> ad.AnnData:
    """Load the Peaks-typed features from a Cell Ranger ARC .h5 as AnnData."""
    a = sc.read_10x_h5(str(path), gex_only=False)
    if "feature_types" not in a.var.columns:
        raise ValueError(f"{path}: feature_types not in .var")
    a = a[:, a.var["feature_types"] == "Peaks"].copy()
    a.var_names_make_unique()
    return a


# ---------------------------------------------------------------------------
# ATAC fragments validation
# ---------------------------------------------------------------------------

def _tabix_list_chromosomes(path: Path) -> list[str] | None:
    """Return chromosome list from the .tbi index via `tabix -l`, or None on failure.

    This is O(1) vs scanning the (potentially huge) fragments file. Requires the
    `tabix` binary on PATH.
    """
    try:
        r = subprocess.run(["tabix", "-l", str(path)], capture_output=True, text=True, check=False)
        if r.returncode != 0:
            return None
        return [c.strip() for c in r.stdout.splitlines() if c.strip()]
    except FileNotFoundError:
        return None


def validate_fragments(path: Path | str, peek_lines: int = 2000) -> dict[str, Any]:
    """Validate structure of fragments.tsv.gz. Returns diagnostics + chromosome set.

    Uses `tabix -l` (O(1)) to enumerate chromosomes; falls back to scanning a small
    peek of the file if tabix is unavailable. Also reads `peek_lines` lines to check
    row structure (5 columns, start<end).
    """
    p = Path(path)
    tbi = Path(str(p) + ".tbi")
    if not p.exists():
        raise FileNotFoundError(p)
    if not tbi.exists():
        raise FileNotFoundError(f"tabix index required: {tbi}")

    # Prefer tabix -l for full chromosome inventory
    tabix_chroms = _tabix_list_chromosomes(p)
    chroms: set[str] = set(tabix_chroms) if tabix_chroms else set()

    # Structural peek — verify first rows look like BED5+
    n_lines = 0
    bad = 0
    with gzip.open(p, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                bad += 1
                continue
            chrom, start, end, _bc, _count = parts[0], parts[1], parts[2], parts[3], parts[4]
            try:
                s, e = int(start), int(end)
                if s >= e:
                    bad += 1
            except ValueError:
                bad += 1
            chroms.add(chrom)
            n_lines += 1
            if n_lines >= peek_lines:
                break

    # Determine prefix convention from canonical chromosomes only (exclude unplaced
    # scaffolds like GL*, JH*, KI* which often lack `chr` prefix in both conventions).
    canonical = [c for c in chroms
                 if c.startswith("chr") or c in {str(i) for i in range(1, 100)}
                 or c in {"X", "Y", "M", "MT"}]
    if canonical:
        has_chr_prefix = all(c.startswith("chr") for c in canonical)
        no_chr_prefix = all(not c.startswith("chr") for c in canonical)
    else:
        has_chr_prefix = False
        no_chr_prefix = True
    if canonical and not (has_chr_prefix or no_chr_prefix):
        raise ValueError(f"Inconsistent chromosome naming in {p}: mix of chr/non-chr on canonical chroms")
    return {
        "path": str(p),
        "tbi_path": str(tbi),
        "peek_lines": n_lines,
        "bad_lines": bad,
        "chromosomes": sorted(chroms),
        "chromosome_source": "tabix" if tabix_chroms else "peek",
        "has_chr_prefix": has_chr_prefix,
    }


GENOME_CHROMS = {
    # minimal fingerprint — presence checks, not exhaustive.
    # We accept BOTH UCSC (chr-prefixed) and Ensembl/NCBI (no-prefix) naming.
    "mm10":   {"chr1", "chr19", "chrX", "chrY", "chrM",
                "1", "19", "X", "Y", "M", "MT"},
    "GRCh38": {"chr1", "chr22", "chrX", "chrY", "chrM",
                "1", "22", "X", "Y", "M", "MT"},
}


def cross_check_genome(chroms: set[str], assembly: str) -> tuple[bool, str]:
    expected = GENOME_CHROMS.get(assembly)
    if expected is None:
        return True, f"No fingerprint for assembly {assembly!r}; skipping cross-check"
    overlap = expected & set(chroms)
    if len(overlap) >= 3:
        return True, f"Found {len(overlap)}/{len(expected)} expected chroms for {assembly}"
    return False, (f"Assembly {assembly} fingerprint mismatch: expected any of "
                   f"{sorted(expected)[:8]}..., got chroms {sorted(chroms)[:8]}...")


# ---------------------------------------------------------------------------
# Fragment barcode extraction (for pairing)
# ---------------------------------------------------------------------------

def fragment_barcodes(path: Path | str, limit: int | None = None) -> set[str]:
    """Scan fragments.tsv.gz and collect the set of barcodes (column 4)."""
    out: set[str] = set()
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                out.add(parts[3])
            if limit is not None and i >= limit:
                break
    return out
