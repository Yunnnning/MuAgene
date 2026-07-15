"""Cleanup / reset policy — the single authority on which working files get deleted, when.

Two phases, both deleting only UNTRACKED working files (never a declared Snakemake
output, so deletion never triggers a re-run):
  - cleanup_qc_intermediates  — at post_qc_review approval (QC/ingest caches).
  - cleanup_process_intermediates — at finish-cleanup, once the S8 deliverable exists.
The durable per-stage markers (validation_report.json / *summary.json / s8_done.txt)
always survive, so `executor status` and the DAG edges hold.
"""
from __future__ import annotations

from pathlib import Path

import yaml

from . import provenance
from .run_paths import RunPaths


# S1a recompute caches: regenerated whenever S1a re-runs, so they are both the
# stale-on-revise set (see revision._S1A_QC_ARTIFACTS) AND safe to delete once QC is
# approved (no further S1a re-run can occur). Single source of truth for both.
S1A_REGEN_CACHES = [
    ("s1a_ambient", "tsne_coords_cache.parquet"),
    ("s1a_ambient", "cell_totals.parquet"),
]


def run_config(run_dir: Path) -> dict:
    """Best-effort load of the run's canonical run.yaml (returns {} on any failure)."""
    try:
        return yaml.safe_load(RunPaths(run_dir).run_yaml.read_text()) or {}
    except Exception:
        return {}


def retain_for_integration(run_dir: Path) -> bool:
    """Whether to KEEP the prepared ATAC fragment caches past the post_qc gate.

    Default True: Integration-MuAgent re-counts fragments against a consensus peak
    set, so the caches must survive approval. A single-sample-and-done run can set
    `retain_for_integration: false` in run.yaml to delete them and reclaim the disk.
    """
    return bool(run_config(run_dir).get("retain_for_integration", True))


def cleanup_qc_intermediates(run_dir: Path) -> list[str]:
    """Remove large QC-only working files after post_qc_review is approved.

    This is the single authority that deletes QC intermediate working files. None
    of these are declared Snakemake outputs — every QC/ingest stage declares only its
    durable marker (s0 validation_report.json, s1a summary.json, s1/s2 qc_summary.json)
    that carries the dependency edges, so removing them here does NOT make any rule's
    declared output "missing" and therefore never triggers a re-run of S0/S1a/S1/S2/S3
    on a later submit. None is read by any post-gate stage (S4–S8); the post-QC h5mu is
    the canonical store S4/S5 read. Removed:
      - rna_qc.h5ad / atac_qc.h5ad — untracked QC matrices, consumed only by S3.
      - atac_snap.h5ad — the pre-filter SnapATAC2 import.
      - atac_snap_explore.h5ad — the qc_explore ATAC import (reused by S2 during
        QC, no longer needed once QC is approved).
      - rna_ingest.h5ad / metadata_minimal.tsv — the S0 RNA ingest (~200 MB) and the
        unused reconstructed metadata TSV. rna_ingest.h5ad is consumed only by S1a (read
        by path) and is a deterministic function of the original input, so it is a pure
        cache: if a re-process needs it again, S1a reconstructs it via io.load_rna_ingest
        (no S0 re-run). metadata_minimal.tsv has no reader. Safe to delete.
      - rna_decontaminated.h5ad — the S1a ambient-corrected RNA (~400 MB), consumed
        only by S1 (read by path). Dead once QC is approved (S4 reads the post-QC h5mu).
      - atac_fragments_cbf[_chrnorm].tsv.gz (+ .tbi) — the chr-normalised fragment
        caches written by io.prepare_fragments_for_snapatac. These are the single
        biggest QC artifact and are RETAINED by default: Integration-MuAgent
        re-counts them against a consensus peak set (reading their *contents*, not
        just the cached filename recorded in parameters.yaml), so the old "dead once
        QC is approved" assumption no longer holds. They are deleted here only when
        run.yaml sets `retain_for_integration: false` (single-sample run reclaiming disk).
      - tsne_coords_cache.parquet / cell_totals.parquet — S1a recompute caches.

    Keeps the durable markers (validation_report.json — also read post-gate by S5 —,
    summary.json, qc_summary.json), the s1/s2 qc_metrics parquets (the durable
    QC-metrics record, consumed by the QC-review summary), and all S3+ artifacts so
    downstream stages are unaffected. Deleting an already-absent file is a no-op.

    Re-processing note: a *post-approval* re-process (re-revise QC thresholds and re-run)
    works without an S0 re-run — S1a reconstructs rna_ingest.h5ad from the original input
    via io.load_rna_ingest, and rna_decontaminated.h5ad regenerates when S1a re-runs. With
    `retain_for_integration` (default true) the fragment caches also survive, so S2 re-runs
    too. Only the previews (figures/plan/qc_explore.json) regenerate from the retained
    *_qc_metrics.parquet (a `revise` does this automatically).
    """
    rp = RunPaths(run_dir)
    # Untracked QC working matrices + the S0/S1a heavy RNA caches (DAG edges are the
    # durable markers validation_report.json / summary.json / qc_summary.json).
    targets = [
        rp.artifact("s1_rna_qc",  "rna_qc.h5ad"),
        rp.artifact("s2_atac_qc", "atac_qc.h5ad"),
        rp.artifact("s2_atac_qc", "atac_snap.h5ad"),
        rp.artifact("qc_explore", "atac_snap_explore.h5ad"),
        rp.artifact("s0_ingest",  "rna_ingest.h5ad"),
        rp.artifact("s0_ingest",  "metadata_minimal.tsv"),
        rp.artifact("s1a_ambient", "rna_decontaminated.h5ad"),
    ]
    # Chr-normalised fragment caches: both naming variants + tabix index, in
    # whichever stage dir they landed in (qc_explore import or an S2 re-derive).
    # RETAINED by default — Integration-MuAgent re-counts them against a consensus
    # peak set, so it reads their contents (not just the cached filename). Deleted
    # only when run.yaml opts out via `retain_for_integration: false`.
    if not retain_for_integration(run_dir):
        for stage in ("qc_explore", "s2_atac_qc"):
            for name in ("atac_fragments_cbf_chrnorm.tsv.gz", "atac_fragments_cbf.tsv.gz"):
                targets.append(rp.artifact(stage, name))
                targets.append(rp.artifact(stage, name + ".tbi"))
    # S1a recompute caches (no S1a re-run can happen after approval).
    targets += [rp.artifact(s, f) for (s, f) in S1A_REGEN_CACHES]

    deleted: list[str] = []
    for p in targets:
        if p.exists():
            p.unlink()
            deleted.append(str(p))
    return deleted


# S4–S8 working files that `finish-cleanup` deletes once the run's final processed
# deliverable exists. All are content-duplicates of the processed h5mu / h5ads
# (normalized RNA, PCA/neighbors, Leiden labels, ATAC spectral embedding, the
# exported peak/tile matrix + its peak-export scratch) — nothing downstream re-reads
# them. None is a declared Snakemake output, and every S4..S8 rule edge depends on the
# durable *_summary.json / s8_done.txt markers (see PROCESS_MARKERS) instead, so
# deleting them never triggers a re-run on a later `submit --target all`.
PROCESS_INTERMEDIATES = [
    # RNA chain (empty stubs on atac_only are removed too)
    ("s4_rna_norm", "rna_norm.h5ad"),
    ("s6_neighbors", "rna_neighbors.h5ad"),
    ("s7_clustering", "rna_clustered.h5ad"),
    # ATAC chain: S5 spectral working file + exported feature sidecars + peak-export
    # scratch (peak_matrix_*.h5ad / *_prepared.bed), and S7 ATAC labels.
    ("s5_atac_spectral", "atac_spectral.h5ad"),
    ("s5_atac_spectral", "feature_matrix.npz"),
    ("s5_atac_spectral", "feature_names.tsv"),
    ("s5_atac_spectral", "feature_kind.txt"),
    ("s5_atac_spectral", "peak_matrix_s2peaks.h5ad"),
    ("s5_atac_spectral", "peak_matrix_user.h5ad"),
    ("s5_atac_spectral", "_s2_peaks_prepared.bed"),
    ("s5_atac_spectral", "_user_peaks_prepared.bed"),
    ("s7_clustering", "atac_leiden_labels.parquet"),
]

# Durable per-stage done-markers that MUST survive finish-cleanup: `executor status`
# and the S4..S8 Snakemake edges key off them. (stage, filename).
PROCESS_MARKERS = [
    ("s4_rna_norm", "norm_summary.json"),
    ("s5_atac_spectral", "spectral_summary.json"),
    ("s6_neighbors", "neighbors_summary.json"),
    ("s7_clustering", "clustering_summary.json"),
    ("s8_umap", "s8_done.txt"),
]


def cleanup_process_intermediates(run_dir: Path) -> list[str]:
    """Remove the large S4–S8 working files once the processed deliverable exists.

    Mirrors `cleanup_qc_intermediates` for the finish phase. Branch-awareness is by
    delete-if-exists: a branch simply never wrote the files it does not apply to
    (rna_only writes no ATAC sidecars), while an atac_only run's empty RNA stubs ARE
    removed. None of these is a declared Snakemake output — the durable markers in
    `PROCESS_MARKERS` carry status + the DAG edges — so deletion never triggers a
    re-run. Deleting an absent file is a no-op. Returns the paths actually deleted.
    """
    rp = RunPaths(run_dir)
    deleted: list[str] = []
    for stage, name in PROCESS_INTERMEDIATES:
        p = rp.artifact(stage, name)
        if p.exists():
            p.unlink()
            deleted.append(str(p))
    return deleted


def processed_outputs_for_branch(rp: RunPaths, branch: str) -> list[Path]:
    """The final S8 processed deliverable(s) for `branch` — what finish-cleanup validates."""
    if branch == "paired":
        return [rp.processed_h5mu]
    out: list[Path] = []
    if branch in ("unpaired", "rna_only"):
        out.append(rp.rna_processed_h5ad)
    if branch in ("unpaired", "atac_only"):
        out.append(rp.atac_processed_h5ad)
    return out


def s8_outputs_valid(run_dir: Path) -> tuple[bool, list[str]]:
    """Whether S8 produced its final processed deliverable(s) — the precondition for
    finish-cleanup. The processed h5mu / h5ads are the user-facing S8 output; a run
    that failed before producing them must keep its intermediates so it can resume.

    Primary check: the branch-derived processed file(s) exist and are non-empty.
    `run_manifest.json` is NOT required (a run can be S8-complete but not yet
    manifest-finalized — e.g. interrupted before the manifest rule); when present, the
    processed paths it records are additionally validated. Returns (ok, problems).
    """
    rp = RunPaths(run_dir)
    branch = provenance.current_branch(str(rp.parameters_yaml))
    problems: list[str] = []

    expected = processed_outputs_for_branch(rp, branch)
    if not expected:
        problems.append(f"unknown workflow_branch {branch!r}; cannot locate S8 outputs")
    for p in expected:
        if not p.exists():
            problems.append(f"missing S8 output: {p}")
        elif p.stat().st_size == 0:
            problems.append(f"empty S8 output: {p}")

    manifest = rp.run_manifest_json
    if manifest.exists():
        import json
        try:
            outputs = (json.loads(manifest.read_text()) or {}).get("outputs") or {}
        except (json.JSONDecodeError, OSError) as e:
            problems.append(f"unreadable run_manifest.json: {e}")
            outputs = {}
        # outputs is a dict of {key: rel_path}; `figures` is a list (non-critical).
        for key, val in outputs.items():
            if key == "figures" or not isinstance(val, str):
                continue
            mp = run_dir / val
            if not mp.exists():
                problems.append(f"manifest output missing: {mp}")
            elif mp.is_file() and mp.stat().st_size == 0:
                problems.append(f"manifest output empty: {mp}")

    return (not problems, problems)


def ensure_process_markers(run_dir: Path) -> list[str]:
    """Backfill any missing S4–S8 durable marker so `executor status` stays `done`
    after cleanup. Safe because the caller only invokes this once the processed
    deliverable is validated (the whole pipeline produced its final output), so every
    stage necessarily completed. Covers legacy / un-finalized runs predating the
    marker refactor (which lack norm/neighbors/clustering_summary.json and/or
    s8_done.txt). Returns the markers written.
    """
    import json
    rp = RunPaths(run_dir)
    written: list[str] = []
    for stage, name in PROCESS_MARKERS:
        p = rp.artifact(stage, name)
        if p.exists():
            continue
        p.parent.mkdir(parents=True, exist_ok=True)
        if name.endswith(".json"):
            p.write_text(json.dumps({"stage": stage, "backfilled": True}, indent=2))
        else:
            p.write_text("backfilled\n")
        written.append(str(p))
    return written
