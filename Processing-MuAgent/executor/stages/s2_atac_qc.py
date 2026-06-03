"""S2 — ATAC QC via SnapATAC2.

Imports fragments into a SnapATAC2-backed AnnData and applies per-cell filters:
  - n_fragments (MAD bounds on log-scale, after an absolute floor)
  - TSS enrichment (min and max)
  - nucleosome signal (max)

Note: tile-matrix construction is NOT part of S2. It happens later in S5 alongside
spectral embedding. S2 only computes QC metrics and subsets cells.
"""
from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from ..methods import mad_thresholds as _mad
from .. import io as _io
from .. import provenance as _prov
from ..log import log_event


def _resolve_param(params_path: Path, plan_params: dict, name: str, default: Any = None) -> Any:
    """parameters.yaml wins over plan (so `executor revise` takes effect on re-run)."""
    v = _prov.get_value(params_path, f"s2_atac_qc.{name}", None)
    if v is not None:
        return v
    entry = plan_params.get(name, {})
    if isinstance(entry, dict) and "value" in entry:
        return entry["value"]
    return default


def _col_to_numpy(adata, key: str) -> np.ndarray:
    """SnapATAC2 obs is polars; convert one column to a numpy float array, or empty array."""
    try:
        col = adata.obs[key]
    except (KeyError, Exception):
        return np.array([], dtype=float)
    # polars Series has .to_numpy(); fall back to np.asarray
    try:
        arr = col.to_numpy()
    except AttributeError:
        arr = np.asarray(col)
    return np.asarray(arr, dtype=float)


def run(run_dir: Path | str, plan: dict[str, Any]) -> dict[str, Any]:
    import snapatac2 as snap
    run_dir = Path(run_dir)
    art = run_dir / "internal" / "artifacts" / "s2_atac_qc"
    art.mkdir(parents=True, exist_ok=True)
    params_path = run_dir / "internal" / "parameters.yaml"

    atac_meta = json.loads((run_dir / "internal" / "artifacts" / "s0_ingest" / "atac_ingest.json").read_text())
    fragments_path = atac_meta["fragments_path"]

    # Optional barcode-translation shim: when S0 persisted a translation parquet
    # (paired branch established via `pairing.translation_table`), stream-rewrite
    # ATAC barcodes into RNA-space before SnapATAC2 ever reads the file. The QC
    # body below is unchanged — it just imports from a different fragments file.
    translation_parquet = (run_dir / "internal" / "artifacts" / "s0_ingest"
                            / "barcode_translation.parquet")
    if translation_parquet.exists():
        from .. import translation as _translation
        translated_path = art / "atac_fragments.translated.tsv.gz"
        if not translated_path.exists():
            table = _translation.load_translation_parquet(translation_parquet)
            stats = _translation.translate_fragments_file(
                fragments_path, translated_path, table,
            )
            log_event(run_dir, {"stage": "s2_atac_qc",
                                 "event": "fragments_translated",
                                 "src": str(fragments_path),
                                 "dst": str(translated_path),
                                 **stats})
        fragments_path = str(translated_path)

    # Strict: require an explicit, SnapATAC2-supported genome. 
    genome = _prov.get_value(params_path, "ingest.genome_assembly", None)
    if not genome:
        raise ValueError(
            "S2 ATAC QC: `ingest.genome_assembly` is not set. Supply it via run.yaml "
            "(genome_assembly: mm10 | GRCh38 | ...) — S2 refuses to guess."
        )
    genome_ref = getattr(snap.genome, genome, None)
    if genome_ref is None:
        available = sorted(n for n in dir(snap.genome) if not n.startswith("_"))
        raise ValueError(
            f"S2 ATAC QC: genome {genome!r} is not supported by SnapATAC2. "
            f"Available assemblies: {available}."
        )

    # Detect Ensembl vs UCSC chromosome naming mismatch between the fragments
    # file and the SnapATAC2 genome object. Cell Ranger ARC (GRCm39 / GRCh38)
    # writes Ensembl-style names ("1", "X") while SnapATAC2 built-in genomes use
    # UCSC-style ("chr1", "chrX"). When the mismatch is detected, add_chr_prefix
    # is passed to filter_fragments_to_chrom_bounds so the output tsv.gz carries
    # UCSC-style names that SnapATAC2 can match.
    genome_uses_chr = any(k.startswith("chr") for k in genome_ref.chrom_sizes)
    frags_chroms = _io._tabix_list_chromosomes(Path(fragments_path))
    add_chr_prefix = False
    if genome_uses_chr and frags_chroms:
        canonical = [c for c in frags_chroms
                     if c in {str(i) for i in range(1, 100)} | {"X", "Y", "M", "MT"}]
        if canonical and not any(c.startswith("chr") for c in canonical):
            add_chr_prefix = True
            log_event(run_dir, {"stage": "s2_atac_qc", "event": "chr_prefix_normalization",
                                 "fragments_chrom_style": "ensembl_no_prefix",
                                 "genome_chrom_style": "ucsc_chr_prefix",
                                 "action": "adding_chr_prefix_in_cbf_filter"})

    # Clip fragments to declared chromosome bounds before SnapATAC2 import.
    # Some aligners produce fragments that extend past a chromosome end
    # (aligner treats end as open interval); SnapATAC2's Rust backend panics on
    # these. Filter is idempotent: skipped if the _cbf file already exists.
    cbf_suffix = "atac_fragments_cbf_chrnorm.tsv.gz" if add_chr_prefix else "atac_fragments_cbf.tsv.gz"
    cbf_path = art / cbf_suffix
    try:
        fragments_path = str(_io.filter_fragments_to_chrom_bounds(
            fragments_path, dict(genome_ref.chrom_sizes), cbf_path,
            add_chr_prefix=add_chr_prefix,
        ))
        _prov.set_param(params_path, "s2_atac_qc.chrom_bound_filter",
                        str(cbf_path),
                        source="derived", confidence="high",
                        rationale=("Fragments with end > chromosome length removed before "
                                   "SnapATAC2 import. Typically <2% of fragments; artifacts of "
                                   "aligners treating chromosome ends as open intervals."),
                        method={"name": "io.filter_fragments_to_chrom_bounds",
                                "code_ref": "executor/io.py::filter_fragments_to_chrom_bounds"})
    except Exception as _e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "chrom_bound_filter_failed",
                             "error": str(_e), "falling_back_to": "unfiltered"})

    # Import fragments into a fresh SnapATAC2-backed h5ad
    h5_out = art / "atac_snap.h5ad"
    adata = snap.pp.import_fragments(
        fragments_path,
        chrom_sizes=genome_ref,
        file=str(h5_out),
        sorted_by_barcode=False,
    )

    # TSS enrichment (per-cell)
    try:
        snap.metrics.tsse(adata, genome_ref)
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "tsse_failed", "error": str(e)})

    n_frag_values = _col_to_numpy(adata, "n_fragment")
    tss_values = _col_to_numpy(adata, "tsse")

    # --- Nucleosome signal (per-cell) ------------------------------------
    # Walk fragments.tsv.gz once, counting nucleosome-free (<147bp) vs
    # mono-nucleosome (147..294bp) fragments per barcode. Signac's NS is
    # mono / nfree. We restrict to S2's import-stage cell set to avoid
    # touching the millions of empty-droplet barcodes.
    cell_barcodes = list(adata.obs_names)
    try:
        ns_values = _io.nucleosome_signal_per_cell(fragments_path, cell_barcodes)
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "nuc_signal_failed",
                            "error": str(e), "fallback": "all_zeros"})
        ns_values = np.zeros(adata.n_obs, dtype=float)

    # Dataset-level fragment-size distribution (cheap; for a sanity figure).
    try:
        snap.metrics.frag_size_distr(adata, max_recorded_size=1000)
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "frag_size_distr_failed",
                            "error": str(e)})

    plan_params = plan["stages"]["s2_atac_qc"]["parameters"]
    k_mad = _resolve_param(params_path, plan_params, "n_fragments_k_mad", 5.0)
    n_frag_floor = _resolve_param(params_path, plan_params, "n_fragments_floor", 500)
    tss_min = float(_resolve_param(params_path, plan_params, "tss_enrichment_min", 1.5))
    tss_max = float(_resolve_param(params_path, plan_params, "tss_enrichment_max", 50.0))
    nuc_signal_max = float(_resolve_param(params_path, plan_params, "nucleosome_signal_max", 3.0))

    if n_frag_values.size:
        keep_floor = n_frag_values >= n_frag_floor
        if keep_floor.any():
            f_lo, f_hi = _mad.log_mad_bounds(n_frag_values[keep_floor], k=k_mad)
        else:
            f_lo, f_hi = float(n_frag_floor), float(n_frag_values.max() if n_frag_values.size else 1e6)
    else:
        f_lo, f_hi = float(n_frag_floor), 1e12

    _prov.set_param(params_path, "s2_atac_qc.n_fragments_min", float(f_lo),
                    source="derived", confidence="high",
                    rationale="MAD lower bound on log1p(n_fragments) after applying floor",
                    method={"name": "mad_thresholds.log_mad_bounds",
                            "code_ref": "executor/methods/mad_thresholds.py"})
    _prov.set_param(params_path, "s2_atac_qc.n_fragments_max", float(f_hi),
                    source="derived", confidence="high",
                    rationale="MAD upper bound on log1p(n_fragments)",
                    method={"name": "mad_thresholds.log_mad_bounds",
                            "code_ref": "executor/methods/mad_thresholds.py"})
    _prov.set_param(params_path, "s2_atac_qc.tss_enrichment_min", float(tss_min),
                    source="recommended", confidence="high",
                    rationale="Minimum TSS enrichment for retained cells")
    _prov.set_param(params_path, "s2_atac_qc.tss_enrichment_max", float(tss_max),
                    source="recommended", confidence="high",
                    rationale="Maximum TSS enrichment; very high values often indicate artifacts")
    _prov.set_param(params_path, "s2_atac_qc.nucleosome_signal_max",
                    float(nuc_signal_max), source="recommended", confidence="high",
                    rationale=("Signac-style NS = mono_nucleosome / nucleosome_free; "
                               "values above this flag poor nucleosome positioning."),
                    method={"name": "io.nucleosome_signal_per_cell",
                            "code_ref": "executor/io.py::nucleosome_signal_per_cell"})

    # Persist per-cell metrics on the AnnData so downstream stages and the
    # qc_summary can read them without re-scanning the fragments file.
    try:
        adata.obs["nucleosome_signal"] = ns_values
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "nuc_signal_obs_write_failed",
                            "error": str(e)})

    keep = np.ones(adata.n_obs, dtype=bool)
    if n_frag_values.size:
        keep &= (n_frag_values >= f_lo) & (n_frag_values <= f_hi)
    if tss_values.size:
        keep &= (tss_values > tss_min) & (tss_values < tss_max)
    if ns_values.size and np.isfinite(ns_values).any():
        keep &= ns_values < nuc_signal_max

    keep_idx = np.nonzero(keep)[0].tolist()
    n_pre = int(adata.n_obs)
    filtered_path = art / "atac_qc.h5ad"
    fd, tmp_name = tempfile.mkstemp(suffix=".h5ad", dir="/tmp")
    Path(tmp_name).unlink(missing_ok=True)
    Path(tmp_name).parent.mkdir(parents=True, exist_ok=True)
    # Close the fd from mkstemp; SnapATAC2 creates the HDF5 file itself.
    import os
    os.close(fd)
    # SnapATAC2 2.8: inplace defaults to True (returns None); pass inplace=False to get
    # a new on-disk AnnData written to `out`. Write on local /tmp first; copying
    # the completed file back to NFS avoids HDF5/NFS locks and Snakemake utime stalls.
    adata_f = adata.subset(obs_indices=keep_idx, out=tmp_name, inplace=False)
    if adata_f is None:
        import snapatac2 as snap
        adata_f = snap.read(tmp_name)
    n_post = int(adata_f.n_obs)

    if n_post == 0:
        try:
            adata_f.close()
        except Exception:
            pass
        Path(tmp_name).unlink(missing_ok=True)
        raise ValueError(
            f"S2 ATAC QC removed all cells (n_pre={n_pre}, n_post=0). Thresholds used: "
            f"n_fragments in [{f_lo:.1f}, {f_hi:.1f}], tss_enrichment in ({tss_min}, {tss_max}), "
            f"nucleosome_signal < {nuc_signal_max}. "
            "Revise one or more thresholds via `executor revise s2_atac_qc ...` "
            "before continuing; downstream stages cannot run on an empty ATAC cell set."
        )

    # Persist summary (+ per-cell ns range for the QC report)
    ns_summary: dict[str, Any] = {}
    if ns_values.size:
        finite_ns = ns_values[np.isfinite(ns_values)]
        if finite_ns.size:
            ns_summary = {
                "median": float(np.median(finite_ns)),
                "p90":    float(np.quantile(finite_ns, 0.90)),
                "max":    float(np.max(finite_ns)),
            }
    _io.write_text_safe(art / "qc_summary.json", json.dumps({
        "n_cells_pre": n_pre,
        "n_cells_post": n_post,
        "thresholds": {
            "n_fragments": [float(f_lo), float(f_hi)],
            "tss_min": float(tss_min),
            "tss_max": float(tss_max),
            "nucleosome_signal_max": float(nuc_signal_max),
        },
        "nucleosome_signal_summary": ns_summary,
    }, indent=2))

    # Dataset-level fragment-size distribution figure — surfaces nucleosome
    # periodicity (well-prepared ATAC libraries show clear ~150 / ~300 / ~450
    # peaks). Cheap to render and very useful for human review.
    try:
        from .. import figures as _fig
        from ..run_paths import RunPaths
        figs_dir = RunPaths(run_dir).deliv_qc_review
        figs_dir.mkdir(parents=True, exist_ok=True)
        # SnapATAC2's `adata.uns` is a polars-backed PyElemCollection, not a
        # plain dict — `.get()` and `in` may not behave like dict semantics.
        # Try-by-key and fall through silently on any access error.
        fsd = None
        try:
            fsd = np.asarray(adata.uns["frag_size_distr"])
        except Exception:
            fsd = None
        if fsd is not None and fsd.size:
            _fig.plot_fragment_size_distribution(
                fsd, out_dir=figs_dir,
                stem="s2_atac_qc_fragment_size_distribution",
                title="ATAC fragment size distribution")
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "frag_size_plot_failed",
                            "error": str(e)})

    try:
        adata.close()
    except Exception:
        pass
    try:
        adata_f.close()
    except Exception:
        pass
    shutil.copy2(tmp_name, filtered_path)
    _io.sync_path(filtered_path)
    Path(tmp_name).unlink(missing_ok=True)

    log_event(run_dir, {"stage": "s2_atac_qc", "event": "done",
                         "n_cells_pre": n_pre, "n_cells_post": n_post,
                         "thresholds": {"n_fragments": [float(f_lo), float(f_hi)],
                                          "tss_min": float(tss_min),
                                          "tss_max": float(tss_max),
                                          "nucleosome_signal_max": float(nuc_signal_max)}})
    return {"n_cells_pre": n_pre, "n_cells_post": n_post}
