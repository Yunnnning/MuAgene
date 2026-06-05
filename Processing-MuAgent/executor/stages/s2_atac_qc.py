"""S2 — ATAC QC via SnapATAC2.

Imports fragments into a SnapATAC2-backed AnnData and applies per-cell filters:
  - n_fragments (MAD bounds on log-scale, after an absolute floor)
  - TSS enrichment (min and max)
  - nucleosome signal (max)
  - FRiP — Fraction of Reads in Peaks (min; only when a peak set is available)

Peak acquisition for FRiP uses the same priority order as S5 feature export:
  0. User-supplied BED (atac_peaks_path in run.yaml)
  1. Cell Ranger ARC pre-called peaks (single_file_multiome) → peak intervals
     extracted from the ARC h5 and written as peaks_arc.bed
  2. MACS3 called on 3-metric-filtered cells → peaks_macs3.bed

When no peak source is available (all tiers fail), FRiP filtering is skipped
and S2 proceeds with the 3-metric-filtered cell set. S5 reuses the BED files
written here rather than calling MACS3 independently.

Note: tile-matrix construction is NOT part of S2. It happens later in S5 alongside
spectral embedding. S2 only computes QC metrics and subsets cells.
"""
from __future__ import annotations

import json
import os
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


def _subset_tss_profile(
    adata,
    obs_indices: list[int],
    genome_ref,
    tmp_dir: Path,
) -> np.ndarray | None:
    """Re-run SnapATAC2 TSS enrichment on a cell subset; return normalized profile."""
    import snapatac2 as snap

    if not obs_indices:
        return None
    fd, tmp = tempfile.mkstemp(suffix=".h5ad", dir=str(tmp_dir))
    Path(tmp).unlink(missing_ok=True)
    os.close(fd)
    sub = adata.subset(obs_indices=obs_indices, out=tmp, inplace=False)
    if sub is None:
        sub = snap.read(tmp)
    try:
        snap.metrics.tsse(sub, genome_ref)
        return np.asarray(sub.uns["TSS_profile"], dtype=float)
    finally:
        try:
            sub.close()
        except Exception:
            pass
        Path(tmp).unlink(missing_ok=True)


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


def _acquire_peaks_for_frip(
    adata_3m,
    run_dir: Path,
    art: Path,
    *,
    genome_ref,
) -> tuple[Path | None, str]:
    """Try to obtain a peak BED file for FRiP computation.

    Returns (peaks_bed_path, source_label) where source_label is one of
    "user_peaks", "arc_h5", "macs3", or "" (empty = no peaks available).
    peaks_bed_path is None when no source succeeded.
    """
    import snapatac2 as snap

    # Priority 0: user-supplied peaks
    try:
        import yaml as _yaml
        from ..run_paths import RunPaths as _RunPaths
        cfg = _yaml.safe_load(_RunPaths(run_dir).run_yaml.read_text()) or {}
        user_peaks = cfg.get("atac_peaks_path")
        if user_peaks:
            p = Path(user_peaks)
            if not p.exists():
                raise RuntimeError(f"atac_peaks_path={user_peaks} not found on disk.")
            return p, "user_peaks"
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "frip_user_peaks_skipped",
                             "reason": str(e)})

    # Priority 1: Cell Ranger ARC pre-called peak intervals
    try:
        import json as _json
        s0_report = _json.loads(
            (run_dir / "internal" / "artifacts" / "s0_ingest" / "validation_report.json").read_text()
        )
        if s0_report.get("single_file_multiome"):
            import yaml as _yaml
            from ..run_paths import RunPaths as _RunPaths
            cfg = _yaml.safe_load(_RunPaths(run_dir).run_yaml.read_text()) or {}
            rna_path = cfg.get("rna_path")
            if rna_path:
                arc_adata = _io.load_atac_from_10x_h5(rna_path)
                arc_peaks_bed = art / "peaks_arc.bed"
                # var_names are genomic intervals, typically "chr1:1000-2000"
                with open(arc_peaks_bed, "w") as fh:
                    for iv in arc_adata.var_names:
                        iv = str(iv)
                        # Handle both "chr1:1000-2000" and "chr1\t1000\t2000" formats
                        if ":" in iv and "-" in iv:
                            chrom, rest = iv.split(":", 1)
                            start, end = rest.split("-", 1)
                        elif "\t" in iv:
                            parts = iv.split("\t")
                            chrom, start, end = parts[0], parts[1], parts[2]
                        else:
                            raise RuntimeError(f"Unrecognised ARC peak interval format: {iv!r}")
                        fh.write(f"{chrom}\t{start}\t{end}\n")
                n_written = len(arc_adata.var_names)
                if n_written == 0:
                    arc_peaks_bed.unlink(missing_ok=True)
                    raise RuntimeError("ARC h5 contained zero ATAC peaks.")
                return arc_peaks_bed, "arc_h5"
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "frip_arc_peaks_skipped",
                             "reason": str(e)})

    # Priority 2: MACS3 called on the 3-metric-filtered cells
    try:
        macs_tempdir = art / "macs3_tmp"
        macs_tempdir.mkdir(parents=True, exist_ok=True)

        # SnapATAC2 2.8: groupby=None → single polars DataFrame;
        # groupby column → dict[str, DataFrame]. Handle both.
        peaks_out = snap.tl.macs3(
            adata_3m, groupby=None, inplace=False,
            qvalue=0.05, shift=-100, extsize=200,
            tempdir=macs_tempdir, n_jobs=1,
        )
        if peaks_out is None:
            raise RuntimeError("snap.tl.macs3 returned no peaks.")

        if isinstance(peaks_out, dict):
            if not peaks_out:
                raise RuntimeError("snap.tl.macs3 returned an empty peak-set dict.")
            merged = snap.tl.merge_peaks(peaks_out, chrom_sizes=genome_ref)
        else:
            merged = peaks_out

        # Normalise column name casing (MACS3 varies between runs)
        cols_lower = [c.lower() for c in merged.columns]
        def _pick(name: str) -> str:
            return merged.columns[cols_lower.index(name)]
        bed_df = merged.select([_pick("chrom"), _pick("start"), _pick("end")])
        if bed_df.height == 0:
            raise RuntimeError("MACS3 produced zero peaks.")

        peaks_macs3_bed = art / "peaks_macs3.bed"
        bed_df.write_csv(str(peaks_macs3_bed), separator="\t", include_header=False)
        return peaks_macs3_bed, "macs3"
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "frip_macs3_skipped",
                             "reason": str(e)})

    return None, ""


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

    # Import fragments into a fresh SnapATAC2-backed h5ad.
    # For Cell Ranger ARC paired multiome, S0 writes a RNA cell-barcode whitelist
    # so we do not import the millions of empty-droplet barcodes in the fragments file.
    whitelist = atac_meta.get("cell_barcode_whitelist")
    h5_out = art / "atac_snap.h5ad"
    adata = snap.pp.import_fragments(
        fragments_path,
        chrom_sizes=genome_ref,
        file=str(h5_out),
        sorted_by_barcode=False,
        whitelist=whitelist,
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
    # Captured here before adata is closed; used later in the figure block.
    fsd_for_fig: np.ndarray | None = None
    try:
        snap.metrics.frag_size_distr(adata, max_recorded_size=1000)
        fsd_for_fig = np.asarray(adata.uns["frag_size_distr"])
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "frag_size_distr_failed",
                            "error": str(e)})

    plan_params = plan["stages"]["s2_atac_qc"]["parameters"]
    k_mad = _resolve_param(params_path, plan_params, "n_fragments_k_mad", 5.0)
    n_frag_floor = _resolve_param(params_path, plan_params, "n_fragments_floor", 1500)
    tss_min = float(_resolve_param(params_path, plan_params, "tss_enrichment_min", 1.5))
    tss_max = float(_resolve_param(params_path, plan_params, "tss_enrichment_max", 50.0))
    nuc_signal_max = float(_resolve_param(params_path, plan_params, "nucleosome_signal_max", 3.0))
    frip_min = float(_resolve_param(params_path, plan_params, "frip_min", 0.2))

    if n_frag_values.size:
        keep_floor = n_frag_values >= n_frag_floor
        if keep_floor.any():
            f_lo, f_hi = _mad.log_mad_bounds(n_frag_values[keep_floor], k=k_mad)
        else:
            f_lo, f_hi = float(n_frag_floor), float(n_frag_values.max() if n_frag_values.size else 1e6)
    else:
        f_lo, f_hi = float(n_frag_floor), 1e12
    f_lo = max(f_lo, float(n_frag_floor))

    _prov.set_param(params_path, "s2_atac_qc.n_fragments_min", float(f_lo),
                    source="derived", confidence="high",
                    rationale=(f"max(MAD lower bound on log1p(n_fragments), "
                               f"n_fragments_floor={n_frag_floor})"),
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
    _prov.set_param(params_path, "s2_atac_qc.frip_min", float(frip_min),
                    source="recommended", confidence="medium",
                    rationale=("Minimum Fraction of Reads in Peaks per cell. "
                               "Set to 0 to disable. Only applied when a peak set is available."))

    # Persist per-cell metrics on the AnnData so downstream stages and the
    # qc_summary can read them without re-scanning the fragments file.
    try:
        adata.obs["nucleosome_signal"] = ns_values
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "nuc_signal_obs_write_failed",
                            "error": str(e)})

    # --- Step A: 3-metric filter (n_fragments, TSS enrichment, nucleosome signal) ---
    keep_3m = np.ones(adata.n_obs, dtype=bool)
    if n_frag_values.size:
        keep_3m &= (n_frag_values >= f_lo) & (n_frag_values <= f_hi)
    if tss_values.size:
        keep_3m &= (tss_values > tss_min) & (tss_values < tss_max)
    if ns_values.size and np.isfinite(ns_values).any():
        keep_3m &= ns_values < nuc_signal_max

    keep_3m_idx = np.nonzero(keep_3m)[0].tolist()
    fail_3m_idx = np.nonzero(~keep_3m)[0].tolist()
    n_pre = int(adata.n_obs)
    n_after_3m = int(len(keep_3m_idx))

    # TSS profile figure: pass = cells meeting all 3-metric thresholds; fail = the rest
    # of the pre-filter import set 
    tss_prof_pass: np.ndarray | None = None
    tss_prof_fail: np.ndarray | None = None
    try:
        tss_prof_pass = _subset_tss_profile(adata, keep_3m_idx, genome_ref, art)
        tss_prof_fail = _subset_tss_profile(adata, fail_3m_idx, genome_ref, art)
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "tss_profile_failed", "error": str(e)})

    # Write 3-metric-filtered AnnData to the artifact dir (NFS has guaranteed space;
    # /tmp on shared compute nodes is small and fills up with large ATAC datasets).
    fd_3m, tmp_3m = tempfile.mkstemp(suffix=".h5ad", dir=str(art))
    Path(tmp_3m).unlink(missing_ok=True)
    os.close(fd_3m)
    adata_3m = adata.subset(obs_indices=keep_3m_idx, out=tmp_3m, inplace=False)
    if adata_3m is None:
        adata_3m = snap.read(tmp_3m)

    try:
        adata.close()
    except Exception:
        pass

    if n_after_3m == 0:
        Path(tmp_3m).unlink(missing_ok=True)
        raise ValueError(
            f"S2 ATAC QC removed all cells after 3-metric filter (n_pre={n_pre}). "
            f"Thresholds: n_fragments in [{f_lo:.1f}, {f_hi:.1f}], "
            f"tss_enrichment in ({tss_min}, {tss_max}), "
            f"nucleosome_signal < {nuc_signal_max}. "
            "Revise one or more thresholds via `executor revise s2_atac_qc ...` "
            "before continuing."
        )

    # --- Step B: Peak acquisition for FRiP ---
    peaks_bed, peak_source = _acquire_peaks_for_frip(
        adata_3m, run_dir, art, genome_ref=genome_ref,
    )

    _prov.set_param(params_path, "s2_atac_qc.peak_source",
                    peak_source if peak_source else "none",
                    source="derived", confidence="high",
                    rationale=("Peak source used for FRiP computation: "
                               "user_peaks | arc_h5 | macs3 | none (FRiP skipped)."),
                    method={"name": "_acquire_peaks_for_frip",
                            "code_ref": "executor/stages/s2_atac_qc.py::_acquire_peaks_for_frip"})

    # --- Step C + D: FRiP computation and filter ---
    frip_values: np.ndarray | None = None
    frip_applied = False

    if peaks_bed is not None and frip_min > 0.0:
        frip_h5 = art / "_frip_tmp.h5ad"
        # SnapATAC2 make_peak_matrix panics on comment lines (e.g. Cell Ranger ARC
        # BED headers starting with '#'). Strip them into a temp file so any peak
        # source — user-supplied, ARC h5 export, or MACS3 — is safe to pass.
        frip_bed_tmp = art / "_peaks_stripped_tmp.bed"
        try:
            frip_h5.unlink(missing_ok=True)
            with open(peaks_bed) as _src, open(frip_bed_tmp, "w") as _dst:
                for _line in _src:
                    if _line.startswith("#"):
                        continue
                    if add_chr_prefix:
                        # Match the chr-prefix added to fragments; BED chrom is
                        # the first whitespace-delimited field.
                        sep = "\t" if "\t" in _line else " "
                        chrom = _line.split(sep, 1)[0]
                        if chrom and not chrom.startswith("chr"):
                            _line = "chr" + _line
                    _dst.write(_line)
            pm_out = snap.pp.make_peak_matrix(
                adata_3m, peak_file=str(frip_bed_tmp), inplace=False, file=str(frip_h5),
            )
            peak_ad = pm_out if pm_out is not None else snap.read(str(frip_h5))

            reads_in_peaks = np.asarray(peak_ad.X[:].sum(axis=1)).ravel()
            n_frag_3m = _col_to_numpy(adata_3m, "n_fragment")
            frip_values = np.where(n_frag_3m > 0, reads_in_peaks / n_frag_3m, 0.0)

            try:
                peak_ad.close()
            except Exception:
                pass
            frip_h5.unlink(missing_ok=True)

            # Persist FRiP in obs so the final AnnData carries it for QC reporting.
            try:
                adata_3m.obs["frip"] = frip_values
            except Exception as e:
                log_event(run_dir, {"stage": "s2_atac_qc", "event": "frip_obs_write_failed",
                                    "error": str(e)})

            frip_applied = True
        except Exception as e:
            log_event(run_dir, {"stage": "s2_atac_qc", "event": "frip_computation_failed",
                                "error": str(e), "fallback": "frip_filter_skipped"})
            frip_values = None
            frip_applied = False
        finally:
            frip_bed_tmp.unlink(missing_ok=True)
    elif peaks_bed is None:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "frip_skipped_no_peaks"})
    # frip_min == 0 means user explicitly disabled FRiP filtering (no log needed)

    # Build the final cell mask on top of the 3m-filtered set
    if frip_applied and frip_values is not None:
        frip_keep = frip_values >= frip_min
        keep_frip_idx = np.nonzero(frip_keep)[0].tolist()
    else:
        keep_frip_idx = list(range(n_after_3m))

    n_post = int(len(keep_frip_idx))

    if n_post == 0:
        try:
            adata_3m.close()
        except Exception:
            pass
        Path(tmp_3m).unlink(missing_ok=True)
        frip_desc = (f", frip >= {frip_min}" if frip_applied else
                     " (FRiP filter not applied — no peak source available)")
        raise ValueError(
            f"S2 ATAC QC removed all cells (n_pre={n_pre}, n_post=0). Thresholds used: "
            f"n_fragments in [{f_lo:.1f}, {f_hi:.1f}], tss_enrichment in ({tss_min}, {tss_max}), "
            f"nucleosome_signal < {nuc_signal_max}{frip_desc}. "
            "Revise one or more thresholds via `executor revise s2_atac_qc ...` "
            "before continuing; downstream stages cannot run on an empty ATAC cell set."
        )

    # Subset to final cell set; write to artifact dir (NFS has space, /tmp can fill).
    fd_f, tmp_final = tempfile.mkstemp(suffix=".h5ad", dir=str(art))
    Path(tmp_final).unlink(missing_ok=True)
    os.close(fd_f)

    if len(keep_frip_idx) == n_after_3m:
        # No FRiP filter applied — rename the 3m file instead of a redundant copy.
        filtered_path_tmp = tmp_3m
        adata_f = adata_3m
    else:
        adata_f = adata_3m.subset(obs_indices=keep_frip_idx, out=tmp_final, inplace=False)
        if adata_f is None:
            adata_f = snap.read(tmp_final)
        try:
            adata_3m.close()
        except Exception:
            pass
        Path(tmp_3m).unlink(missing_ok=True)
        filtered_path_tmp = tmp_final

    # Persist summary (+ per-cell ns range and FRiP info for the QC report)
    ns_summary: dict[str, Any] = {}
    if ns_values.size:
        finite_ns = ns_values[np.isfinite(ns_values)]
        if finite_ns.size:
            ns_summary = {
                "median": float(np.median(finite_ns)),
                "p90":    float(np.quantile(finite_ns, 0.90)),
                "max":    float(np.max(finite_ns)),
            }

    frip_summary: dict[str, Any] | None = None
    if frip_applied and frip_values is not None and frip_values.size:
        frip_summary = {
            "median": float(np.median(frip_values)),
            "p10":    float(np.quantile(frip_values, 0.10)),
            "min":    float(np.min(frip_values)),
        }

    _io.write_text_safe(art / "qc_summary.json", json.dumps({
        "n_cells_pre": n_pre,
        "n_cells_after_3m_filter": n_after_3m,
        "n_cells_post": n_post,
        "thresholds": {
            "n_fragments": [float(f_lo), float(f_hi)],
            "tss_min": float(tss_min),
            "tss_max": float(tss_max),
            "nucleosome_signal_max": float(nuc_signal_max),
            "frip_min": float(frip_min) if frip_applied else None,
        },
        "nucleosome_signal_summary": ns_summary,
        "frip_summary": frip_summary,
        "peak_source": peak_source if peak_source else None,
    }, indent=2))

    # Capture post-QC fragment size distribution before adata_f is closed.
    fsd_after: np.ndarray | None = None
    try:
        snap.metrics.frag_size_distr(adata_f, max_recorded_size=1000)
        fsd_after = np.asarray(adata_f.uns["frag_size_distr"])
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "frag_size_distr_post_failed",
                             "error": str(e)})

    # Figures: (1) fragment-size distribution, (2) TSS enrichment profile,
    #          (3) FRiP histogram with threshold line.
    try:
        from .. import figures as _fig
        from ..run_paths import RunPaths
        figs_dir = RunPaths(run_dir).deliv_qc_review
        figs_dir.mkdir(parents=True, exist_ok=True)

        if fsd_for_fig is not None and fsd_for_fig.size:
            _fig.plot_fragment_size_distribution(
                fsd_for_fig, out_dir=figs_dir,
                stem="s2_atac_qc_fragment_size_distribution",
                title="ATAC fragment size distribution (after QC)",
                distr_after=fsd_after if (fsd_after is not None and fsd_after.size) else None,
            )

        if tss_prof_pass is not None and tss_prof_fail is not None:
            _fig.plot_tss_enrichment_profile(
                tss_prof_pass, tss_prof_fail,
                out_dir=figs_dir,
                stem="s2_atac_qc_tss_enrichment_profile",
                n_pass=n_after_3m,
                n_fail=n_pre - n_after_3m,
            )

        if frip_applied and frip_values is not None and frip_values.size:
            _fig.plot_frip_histogram(
                frip_values, frip_min=frip_min,
                out_dir=figs_dir,
                stem="s2_atac_qc_frip_histogram",
            )
    except Exception as e:
        log_event(run_dir, {"stage": "s2_atac_qc", "event": "figures_failed",
                             "error": str(e)})

    try:
        adata_f.close()
    except Exception:
        pass

    final_artifact = art / "atac_qc.h5ad"
    shutil.copy2(filtered_path_tmp, final_artifact)
    _io.sync_path(final_artifact)
    Path(filtered_path_tmp).unlink(missing_ok=True)

    log_event(run_dir, {"stage": "s2_atac_qc", "event": "done",
                         "n_cells_pre": n_pre,
                         "n_cells_after_3m_filter": n_after_3m,
                         "n_cells_post": n_post,
                         "peak_source": peak_source if peak_source else None,
                         "frip_applied": frip_applied,
                         "thresholds": {"n_fragments": [float(f_lo), float(f_hi)],
                                          "tss_min": float(tss_min),
                                          "tss_max": float(tss_max),
                                          "nucleosome_signal_max": float(nuc_signal_max),
                                          "frip_min": float(frip_min) if frip_applied else None}})
    return {"n_cells_pre": n_pre, "n_cells_after_3m_filter": n_after_3m, "n_cells_post": n_post}
