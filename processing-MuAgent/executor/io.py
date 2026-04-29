"""RNA/ATAC input format autodetection and loading.

Supported RNA formats: 10x HDF5 (.h5), 10x MEX dir, AnnData (.h5ad), custom.
ATAC: fragments.tsv.gz (+ .tbi required), optional peak matrix / h5ad.

Raw-vs-filtered status is detected via barcode count: 10x raw matrices contain
the full whitelist (~6.7M barcodes); filtered/cell-called matrices contain
~50K cells or fewer. The threshold below is conservative.
"""
from __future__ import annotations

import gzip
import json
import subprocess
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import scanpy as sc


# Above this barcode count the matrix is treated as a raw (cell-not-called)
# matrix. 10x raw outputs typically have ~6.7M barcodes; filtered have <100K.
RAW_BARCODE_THRESHOLD = 200_000


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
# Raw vs filtered RNA matrix detection + barcode-rank cell calling
# ---------------------------------------------------------------------------

def detect_filtered_status(path: Path | str, fmt: str | None = None) -> str:
    """Return "filtered" | "raw" by peeking at the barcode count.

    Threshold is RAW_BARCODE_THRESHOLD; intentionally conservative so a 10x
    raw matrix (~6.7M barcodes) is unambiguously classified.
    """
    p = Path(path)
    fmt = fmt or detect_rna_format(p)
    n_barcodes: int | None = None
    if fmt == "10x_h5":
        try:
            with h5py.File(p, "r") as f:
                if "matrix" in f and "barcodes" in f["matrix"]:
                    n_barcodes = int(f["matrix"]["barcodes"].shape[0])
        except Exception:
            n_barcodes = None
    elif fmt == "10x_mex":
        bc = p / "barcodes.tsv.gz"
        if not bc.exists():
            bc = p / "barcodes.tsv"
        if bc.exists():
            try:
                opener = gzip.open if str(bc).endswith(".gz") else open
                with opener(bc, "rt") as f:
                    n_barcodes = sum(1 for _ in f)
            except Exception:
                n_barcodes = None
    elif fmt == "h5ad":
        try:
            with h5py.File(p, "r") as f:
                if "obs" in f:
                    obs_grp = f["obs"]
                    if "_index" in obs_grp:
                        n_barcodes = int(obs_grp["_index"].shape[0])
                    elif "index" in obs_grp:
                        n_barcodes = int(obs_grp["index"].shape[0])
        except Exception:
            n_barcodes = None
    if n_barcodes is None:
        return "filtered"
    return "raw" if n_barcodes >= RAW_BARCODE_THRESHOLD else "filtered"


def barcode_rank_knee(total_counts: np.ndarray) -> tuple[int, dict[str, Any]]:
    """Knee-point cell calling on the barcode-rank (log-counts vs log-rank) curve.

    Implements a curvature-based knee finder roughly equivalent to the 10x
    barcode-rank "knee" call (and the kneedle algorithm idea), without the
    `kneed` dependency. Returns the count threshold and a small diagnostic
    dict; cells with `total_counts >= threshold` are kept.

    The curve is unimodal-decreasing in log-log; the knee is the rank where
    log(counts) drops most steeply. We pick the rank that maximises the
    distance from each point to the chord between (rank=1, max_count) and
    (rank=N, min_count). Ranks below the knee = real cells.
    """
    counts = np.asarray(total_counts, dtype=float)
    counts = counts[counts > 0]
    if counts.size < 100:
        thresh = float(np.percentile(counts, 50)) if counts.size else 0.0
        return int(np.sum(counts >= thresh)), {"method": "fallback_p50",
                                                "threshold": thresh,
                                                "n_kept": int(np.sum(counts >= thresh))}
    sorted_counts = np.sort(counts)[::-1]
    log_rank = np.log10(np.arange(1, sorted_counts.size + 1))
    log_counts = np.log10(np.maximum(sorted_counts, 1.0))
    p1 = np.array([log_rank[0], log_counts[0]])
    p2 = np.array([log_rank[-1], log_counts[-1]])
    chord = p2 - p1
    chord_norm = np.linalg.norm(chord)
    if chord_norm == 0:
        thresh = float(np.percentile(counts, 50))
        return int(np.sum(counts >= thresh)), {"method": "fallback_degenerate_chord",
                                                "threshold": thresh,
                                                "n_kept": int(np.sum(counts >= thresh))}
    chord_unit = chord / chord_norm
    points = np.column_stack([log_rank, log_counts])
    rel = points - p1
    proj = rel - np.outer(rel @ chord_unit, chord_unit)
    distances = np.linalg.norm(proj, axis=1)
    # Restrict knee search to the upper region (avoid picking the long tail).
    # Use distance only where log_counts > median; fall back if empty.
    upper_mask = log_counts > np.median(log_counts)
    if upper_mask.any():
        masked = np.where(upper_mask, distances, -np.inf)
        knee_idx = int(np.argmax(masked))
    else:
        knee_idx = int(np.argmax(distances))
    threshold = float(sorted_counts[knee_idx])
    n_kept = int(np.sum(counts >= threshold))
    return n_kept, {
        "method": "barcode_rank_knee_chord_distance",
        "threshold": threshold,
        "knee_rank": int(knee_idx + 1),
        "n_kept": n_kept,
        "n_barcodes": int(counts.size),
    }


def call_cells_from_raw(adata: ad.AnnData, *, min_counts_floor: int = 100) -> tuple[ad.AnnData, dict[str, Any]]:
    """Apply barcode-rank knee cell calling to a raw RNA AnnData.

    Returns `(adata_filtered, diag_dict)` where `adata_filtered` contains only
    barcodes with `total_counts >= knee_threshold`. An absolute floor of
    `min_counts_floor` guards against degenerate inputs. Original raw matrix
    is preserved on the caller side; this function does not mutate `adata`
    in-place beyond returning a view-derived copy.
    """
    X = adata.X
    if hasattr(X, "sum"):
        total = np.asarray(X.sum(axis=1)).ravel()
    else:
        total = np.asarray(X).sum(axis=1)
    n_kept, diag = barcode_rank_knee(total)
    threshold = max(diag.get("threshold", 0.0), float(min_counts_floor))
    keep = total >= threshold
    diag["threshold"] = float(threshold)
    diag["n_kept"] = int(keep.sum())
    diag["n_dropped"] = int((~keep).sum())
    diag["min_counts_floor"] = int(min_counts_floor)
    return adata[keep].copy(), diag


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


def fragment_size_metrics(
    path: Path | str,
    barcodes: list[str] | set[str] | None = None,
    *,
    nfree_max: int = 147,
    mono_max: int = 294,
) -> dict[str, dict[str, int]]:
    """Per-cell nucleosome-band fragment counts (one pass over fragments.tsv.gz).

    Bins each fragment by length (`end - start`) into:
      - nucleosome_free:  1 <= L < `nfree_max` (default 147)
      - mono_nucleosome:  `nfree_max` <= L < `mono_max` (default 294)
    Returns `{barcode: {"nfree": int, "mono": int}}`.

    The Signac `nucleosome_signal` is then `mono / nfree` per cell, which is
    what S2 filters on. Fragments outside the two bins (>= mono_max) are
    discarded — they are uninformative for the nucleosome ratio.

    Restricting to a barcode set (typically S2's import-stage cell set) avoids
    materialising counts for the millions of low-quality droplets in raw
    fragments files.
    """
    bc_filter: set[str] | None = None
    if barcodes is not None:
        bc_filter = set(barcodes)
    out: dict[str, dict[str, int]] = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            bc = parts[3]
            if bc_filter is not None and bc not in bc_filter:
                continue
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            length = end - start
            if length < 1 or length >= mono_max:
                continue
            rec = out.setdefault(bc, {"nfree": 0, "mono": 0})
            if length < nfree_max:
                rec["nfree"] += 1
            else:
                rec["mono"] += 1
    return out


def nucleosome_signal_per_cell(
    path: Path | str,
    barcodes: list[str],
    *,
    pseudocount: float = 1.0,
) -> np.ndarray:
    """Compute Signac-style `nucleosome_signal = mono / nfree` per cell.

    `barcodes` defines the order of the returned 1D array; cells with no
    fragments in either band get value `mono / pseudocount` (which is 0.0 by
    default if `mono == 0`). The pseudocount avoids divide-by-zero on cells
    with no nucleosome-free fragments.
    """
    metrics = fragment_size_metrics(path, barcodes=barcodes)
    out = np.zeros(len(barcodes), dtype=float)
    for i, bc in enumerate(barcodes):
        rec = metrics.get(bc, {"nfree": 0, "mono": 0})
        denom = rec["nfree"] + pseudocount
        out[i] = rec["mono"] / denom
    return out
