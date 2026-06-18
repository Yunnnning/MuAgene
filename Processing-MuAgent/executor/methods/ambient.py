"""Ambient-RNA correction wrappers.

Two backends are supported and dispatched by S1a based on what S0 produced:

  - **DecontX** (celda, R) — uses the filtered counts matrix only. Picked when
    no companion raw matrix exists.
  - **SoupX** (R) — uses both the filtered cells and the raw droplets to
    estimate the soup profile. Picked when S0 registers a raw-matrix input
    ref (symlink to `rna_raw_path`, or to `rna_path` when that matrix is raw).

Both are invoked via `Rscript`. The Python side serialises counts to Matrix
Market (.mtx) for portability — neither AnnData round-tripping nor rpy2 is
required at runtime. R, celda, and SoupX are listed in
`workflow/envs/processing.yaml`. If they are missing when correction is
requested, wrappers raise `AmbientUnavailable` and S1a fails the run.
Pass-through is only used when the plan sets `method=none` or the RNA matrix
is empty (e.g. `atac_only` branch).
"""
from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy.io import mmread, mmwrite


R_SCRIPTS_DIR = Path(__file__).resolve().parent / "r_scripts"


class AmbientUnavailable(RuntimeError):
    """R / DecontX / SoupX missing — recreate the conda env from processing.yaml."""


@dataclass
class AmbientResult:
    """Output of an ambient-correction call.

    `corrected_counts` is the cells x genes integer-valued sparse counts matrix
    (cells in the same order as `barcodes`, genes as `features`); `contamination`
    is per-cell rho/contamination fraction in [0, 1]; `summary` is the JSON
    payload written by the R wrapper.
    """

    method: str
    corrected_counts: sp.csr_matrix
    contamination: np.ndarray
    barcodes: list[str]
    features: list[str]
    summary: dict[str, Any]


# ---------------------------------------------------------------------------
# R availability + script invocation
# ---------------------------------------------------------------------------

def rscript_path() -> str | None:
    return shutil.which("Rscript")


def check_r_packages(packages: list[str]) -> dict[str, bool]:
    """Return `{package: installed_bool}` for each package via a single Rscript call.

    Returns all-False when Rscript itself is missing.
    """
    out: dict[str, bool] = {p: False for p in packages}
    if rscript_path() is None:
        return out
    expr = ";".join(
        f'cat("{p}=", as.integer(suppressWarnings(suppressMessages(requireNamespace("{p}", quietly=TRUE)))), "\\n", sep="")'
        for p in packages
    )
    try:
        r = subprocess.run([rscript_path(), "-e", expr],
                            capture_output=True, text=True, check=False, timeout=60)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return out
    for line in r.stdout.splitlines():
        if "=" in line:
            name, val = line.split("=", 1)
            out[name.strip()] = val.strip().startswith("1")
    return out


def _run_rscript(script: Path, args: list[str], *, log_file: Path | None = None) -> None:
    rscript = rscript_path()
    if rscript is None:
        raise AmbientUnavailable(
            "Rscript is not on PATH; activate the muagene env "
            "(workflow/envs/processing.yaml)."
        )
    cmd = [rscript, "--vanilla", str(script), *args]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        log_file.write_text(
            f"$ {' '.join(cmd)}\n"
            f"--- stdout ---\n{proc.stdout}\n"
            f"--- stderr ---\n{proc.stderr}\n"
            f"--- exit_code: {proc.returncode} ---\n"
        )
    if proc.returncode != 0:
        msg = proc.stderr or proc.stdout
        raise AmbientUnavailable(
            f"Rscript {script.name} failed (exit={proc.returncode}). "
            f"See log{f' {log_file}' if log_file else ''}.\n{msg.strip()[-1500:]}"
        )


# ---------------------------------------------------------------------------
# Counts (de)serialisation helpers
# ---------------------------------------------------------------------------

def _to_csr(adata: ad.AnnData) -> sp.csr_matrix:
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)


def _write_mtx(matrix: sp.spmatrix, path: Path) -> None:
    mmwrite(str(path), matrix.tocoo())


def _read_mtx(path: Path) -> sp.csr_matrix:
    return mmread(str(path)).tocsr()


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------

def run_decontx(
    adata: ad.AnnData,
    *,
    work_dir: Path,
    max_contamination: float = 1.0,
) -> AmbientResult:
    """Run DecontX on a filtered counts matrix. Cells in `adata.obs_names`."""
    pkgs = check_r_packages(["celda", "Matrix", "jsonlite"])
    if not all(pkgs.values()):
        missing = [p for p, ok in pkgs.items() if not ok]
        raise AmbientUnavailable(
            f"DecontX requires R packages {missing}; recreate the conda env from "
            "workflow/envs/processing.yaml (bioconductor-celda, r-base)."
        )

    in_dir = work_dir / "in"
    out_dir = work_dir / "out"
    in_dir.mkdir(parents=True, exist_ok=True)

    counts = _to_csr(adata)
    if counts.nnz > 0:
        if counts.data.min() < 0:
            raise ValueError(
                "DecontX requires non-negative counts; negative values found in the count matrix."
            )
        if not np.allclose(counts.data, np.round(counts.data)):
            raise ValueError(
                "DecontX requires integer counts; non-integer values found "
                "(data may be normalized rather than raw counts)."
            )
    barcodes = list(map(str, adata.obs_names))
    features = list(map(str, adata.var_names))
    _write_mtx(counts, in_dir / "counts.mtx")
    (in_dir / "barcodes.tsv").write_text("\n".join(barcodes) + "\n")
    (in_dir / "features.tsv").write_text("\n".join(features) + "\n")

    _run_rscript(
        R_SCRIPTS_DIR / "decontx.R",
        [str(in_dir), str(out_dir), str(max_contamination)],
        log_file=work_dir / "rscript.log",
    )

    decon = _read_mtx(out_dir / "decontaminated.mtx")
    contam = np.loadtxt(out_dir / "contamination.tsv", dtype=float).reshape(-1)
    summary = json.loads((out_dir / "summary.json").read_text())
    return AmbientResult(
        method="DecontX",
        corrected_counts=decon,
        contamination=contam,
        barcodes=barcodes,
        features=features,
        summary=summary,
    )


def run_soupx(
    filtered: ad.AnnData,
    raw: ad.AnnData,
    *,
    work_dir: Path,
    max_contamination: float = 1.0,
) -> AmbientResult:
    """Run SoupX given filtered cells + raw droplets.

    The raw matrix's gene set is intersected with `filtered.var_names`, and the
    filtered matrix is transposed onto that shared gene order before calling R.
    """
    pkgs = check_r_packages(["SoupX", "Matrix", "jsonlite"])
    if not all(pkgs.values()):
        missing = [p for p, ok in pkgs.items() if not ok]
        raise AmbientUnavailable(
            f"SoupX requires R packages {missing}; recreate the conda env from "
            "workflow/envs/processing.yaml (r-soupx, r-base)."
        )

    # Align gene order between raw and filtered (intersection in filtered's order).
    shared = [g for g in filtered.var_names if g in set(raw.var_names)]
    if len(shared) < 200:
        raise ValueError(
            f"SoupX: only {len(shared)} genes shared between raw and filtered "
            "matrices; the matrices are not from the same reference."
        )
    filt_aligned = filtered[:, shared].copy()
    raw_aligned = raw[:, shared].copy()

    in_dir = work_dir / "in"
    out_dir = work_dir / "out"
    in_dir.mkdir(parents=True, exist_ok=True)

    filt_X = _to_csr(filt_aligned)
    raw_X = _to_csr(raw_aligned)
    filt_bc = list(map(str, filt_aligned.obs_names))
    raw_bc = list(map(str, raw_aligned.obs_names))
    features = list(map(str, shared))

    _write_mtx(filt_X, in_dir / "filtered_counts.mtx")
    _write_mtx(raw_X, in_dir / "raw_counts.mtx")
    (in_dir / "filtered_barcodes.tsv").write_text("\n".join(filt_bc) + "\n")
    (in_dir / "raw_barcodes.tsv").write_text("\n".join(raw_bc) + "\n")
    (in_dir / "features.tsv").write_text("\n".join(features) + "\n")

    _run_rscript(
        R_SCRIPTS_DIR / "soupx.R",
        [str(in_dir), str(out_dir), str(max_contamination)],
        log_file=work_dir / "rscript.log",
    )

    decon = _read_mtx(out_dir / "decontaminated.mtx")
    contam = np.loadtxt(out_dir / "contamination.tsv", dtype=float).reshape(-1)
    summary = json.loads((out_dir / "summary.json").read_text())
    return AmbientResult(
        method="SoupX",
        corrected_counts=decon,
        contamination=contam,
        barcodes=filt_bc,
        features=features,
        summary=summary,
    )


def apply_correction(
    filtered: ad.AnnData,
    result: AmbientResult,
) -> ad.AnnData:
    """Return a copy of `filtered` with corrected counts in `.X` + `.layers`.

    Adds:
        - `.layers['counts_raw']`: the original (pre-correction) counts.
        - `.layers['counts']`:     the decontaminated integer counts (also `.X`).
        - `.obs['ambient_contamination']`: per-cell contamination fraction.
        - `.uns['ambient']`:           full method/summary diagnostic dict.

    The cell order is preserved; genes are restricted to `result.features` if
    they're a subset (SoupX path), otherwise all genes are kept (DecontX path).
    """
    out = filtered.copy()
    if list(out.var_names) != list(result.features):
        out = out[:, result.features].copy()
    if list(out.obs_names) != list(result.barcodes):
        out = out[result.barcodes, :].copy()

    raw_layer = out.layers.get("counts", out.X)
    out.layers["counts_raw"] = raw_layer.copy() if hasattr(raw_layer, "copy") else raw_layer
    out.layers["counts"] = result.corrected_counts.tocsr()
    out.X = result.corrected_counts.tocsr()
    out.obs["ambient_contamination"] = np.asarray(result.contamination, dtype=float)
    out.uns["ambient"] = {"method": result.method, **result.summary}
    return out
