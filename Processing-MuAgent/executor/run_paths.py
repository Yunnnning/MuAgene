"""Single source of truth for canonical run-directory paths.

Top-level layout (direct-write; no symlinks):

    <run_dir>/
      deliverables/                 user-facing; split by lifecycle phase
        pre_run/                    materials produced BEFORE preprocessing execution
          config/
            run.yaml                    ← canonical copy of the user config
            biological_context.md       ← canonical Biological Context Report
          summary/
            context_summary.md          ← P1 output
            plan_review.md              ← plan review gate (summary + parameter appendix)
        checkpoint/                 intermediate review artifacts (flat subfolders)
          qc_review/                  QC review checkpoint (figures + qc_review.md)
            qc_review.md
            s1_rna_qc_violin_{pre,post}.{png,pdf}
            ...
          resolution_review/          S7 resolution checkpoint
            resolution_summary.md
            resolution_review.{ipynb,html}
            s7_compare_*.{png,pdf}    (from resolution-compare CLI)
        post_run/                   final deliverables only (flat — no subfolders)
          s8_umap_{rna,atac}_by_leiden.{png,pdf}
          review_processed_h5mu.{ipynb,py}
          processed.h5mu              (paired branch)
          rna_processed.h5ad          (separate / rna_only branches)
          atac_processed.h5ad         (separate / atac_only branches)
          qc_summary.md               ← manifest rule
          run_manifest.json           ← manifest rule (handoff artifact)
          layout.json                 ← layout.finalize (manifest of deliverables)
      internal/                     canonical pipeline state (not user-facing)
        artifacts/sN_<stage>/       intermediate stage outputs
        proposals/                  <stage>.yaml + awaiting_approval sentinels
        checkpoints/                <stage>.approved sentinels
        parameters.yaml
        state.yaml
        log.jsonl

Key invariants:
- No file lives at the top of <run_dir>/; exactly two subdirectories exist there.
- No symlinks or aliases — every logical artifact has a single canonical path.
- QC figures + checkpoint summaries live under checkpoint/{qc_review,resolution_review}/.
- post_run/ contains UMAP figures, processed data, final notebook, and manifest only.
- pre_run/ contains only files the user reviews before approving the plan.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RunPaths:
    """Frozen wrapper around a `run_dir` exposing every canonical path."""

    run_dir: Path

    def __post_init__(self) -> None:
        object.__setattr__(self, "run_dir", Path(self.run_dir))

    # --- Top-level ---------------------------------------------------------
    @property
    def internal(self) -> Path:
        return self.run_dir / "internal"

    @property
    def deliverables(self) -> Path:
        return self.run_dir / "deliverables"

    # --- Internal state (non-user-facing) ---------------------------------
    @property
    def artifacts(self) -> Path:
        return self.internal / "artifacts"

    @property
    def proposals(self) -> Path:
        return self.internal / "proposals"

    @property
    def checkpoints(self) -> Path:
        return self.internal / "checkpoints"

    @property
    def parameters_yaml(self) -> Path:
        return self.internal / "parameters.yaml"

    @property
    def state_yaml(self) -> Path:
        return self.internal / "state.yaml"

    @property
    def log_jsonl(self) -> Path:
        return self.internal / "log.jsonl"

    @property
    def snakemake_workdir(self) -> Path:
        """Run-local Snakemake state directory.

        Snakemake writes locks, metadata, and scheduler logs under `.snakemake`
        in its working directory. Keeping that directory inside each run avoids
        cross-run lock contention when multiple datasets run concurrently.
        """
        return self.internal / "snakemake"

    # --- Deliverable phase roots ------------------------------------------
    @property
    def deliv_pre_run(self) -> Path:
        return self.deliverables / "pre_run"

    @property
    def deliv_checkpoint(self) -> Path:
        return self.deliverables / "checkpoint"

    @property
    def deliv_qc_review(self) -> Path:
        """checkpoint/qc_review/ — QC figures + pre-dimred summary."""
        return self.deliv_checkpoint / "qc_review"

    @property
    def deliv_resolution_review(self) -> Path:
        """checkpoint/resolution_review/ — S7 summary, notebook, comparison figures."""
        return self.deliv_checkpoint / "resolution_review"

    @property
    def deliv_post_run(self) -> Path:
        """post_run/ — flat final deliverables (UMAP figures, data, notebook, manifest)."""
        return self.deliverables / "post_run"

    # --- pre_run sub-directories (unchanged) ------------------------------
    @property
    def deliv_config(self) -> Path:
        """pre_run/config/ — user-supplied config + biological context."""
        return self.deliv_pre_run / "config"

    @property
    def deliv_pre_summary(self) -> Path:
        """pre_run/summary/ — P1 / P2 / plan_review materials for approval."""
        return self.deliv_pre_run / "summary"

    # --- Canonical user-facing files --------------------------------------
    @property
    def run_yaml(self) -> Path:
        return self.deliv_config / "run.yaml"

    @property
    def biological_context_md(self) -> Path:
        return self.deliv_config / "biological_context.md"

    @property
    def hpc_env_sh(self) -> Path:
        """pre_run/config/hpc.env — source-able PMA_* exports for cluster runs."""
        return self.deliv_config / "hpc.env"

    @property
    def site_config(self) -> Path:
        """pre_run/config/site.config — YAML platform description consumed by Execution-MuAgent."""
        return self.deliv_config / "site.config"

    @property
    def stage_meta_dir(self) -> Path:
        """internal/stage_meta/ — monitoring metadata (resources, I/O, timeout hints). Not a submission contract."""
        return self.internal / "stage_meta"

    def stage_meta(self, stage: str) -> Path:
        """internal/stage_meta/<stage>.yaml"""
        return self.stage_meta_dir / f"{stage}.yaml"

    # pre_run/summary/*
    @property
    def context_summary_md(self) -> Path:
        return self.deliv_pre_summary / "context_summary.md"

    @property
    def plan_review_md(self) -> Path:
        """pre_run/summary/plan_review.md — merged summary + full parameter appendix."""
        return self.deliv_pre_summary / "plan_review.md"

    # checkpoint/qc_review/*
    @property
    def qc_review_summary_md(self) -> Path:
        """checkpoint/qc_review/qc_review.md — QC review user checkpoint (after S3)."""
        return self.deliv_qc_review / "qc_review.md"

    @property
    def qc_summary_pre_dimred_md(self) -> Path:
        """Deprecated alias for qc_review_summary_md."""
        return self.qc_review_summary_md

    # checkpoint/resolution_review/*
    @property
    def resolution_summary_md(self) -> Path:
        return self.deliv_resolution_review / "resolution_summary.md"

    @property
    def resolution_review_ipynb(self) -> Path:
        return self.deliv_resolution_review / "resolution_review.ipynb"

    @property
    def resolution_review_html(self) -> Path:
        return self.deliv_resolution_review / "resolution_review.html"

    # post_run/* (flat)
    @property
    def qc_summary_md(self) -> Path:
        return self.deliv_post_run / "qc_summary.md"

    @property
    def run_manifest_json(self) -> Path:
        return self.deliv_post_run / "run_manifest.json"

    @property
    def layout_json(self) -> Path:
        return self.deliv_post_run / "layout.json"

    @property
    def processed_h5mu(self) -> Path:
        return self.deliv_post_run / "processed.h5mu"

    @property
    def rna_processed_h5ad(self) -> Path:
        return self.deliv_post_run / "rna_processed.h5ad"

    @property
    def atac_processed_h5ad(self) -> Path:
        return self.deliv_post_run / "atac_processed.h5ad"

    @property
    def review_notebook_ipynb(self) -> Path:
        return self.deliv_post_run / "review_processed_h5mu.ipynb"

    @property
    def review_notebook_py(self) -> Path:
        return self.deliv_post_run / "review_processed_h5mu.py"

    # --- Stage helpers -----------------------------------------------------
    def stage_dir(self, stage: str) -> Path:
        """Return `internal/artifacts/<stage>/`."""
        return self.artifacts / stage

    def artifact(self, stage: str, *parts: str) -> Path:
        """Return `internal/artifacts/<stage>/<parts>`."""
        return self.stage_dir(stage).joinpath(*parts)

    def proposal(self, stage: str, suffix: str = ".yaml") -> Path:
        return self.proposals / f"{stage}{suffix}"

    def awaiting_sentinel(self, stage: str) -> Path:
        return self.proposals / f"{stage}.awaiting_approval"

    def approved_sentinel(self, stage: str) -> Path:
        return self.checkpoints / f"{stage}.approved"

    def deliv_qc_figure(self, stem: str, *, ext: str = "png") -> Path:
        return self.deliv_qc_review / f"{stem}.{ext}"

    def deliv_umap_figure(self, stem: str, *, ext: str = "png") -> Path:
        return self.deliv_post_run / f"{stem}.{ext}"

    # --- Bootstrap ---------------------------------------------------------
    def ensure(self) -> None:
        """Create both the internal and deliverables scaffold (idempotent)."""
        for p in (
            self.internal,
            self.artifacts, self.proposals, self.checkpoints, self.snakemake_workdir,
            self.deliverables,
            self.deliv_pre_run, self.deliv_config, self.deliv_pre_summary,
            self.deliv_checkpoint, self.deliv_qc_review, self.deliv_resolution_review,
            self.deliv_post_run,
        ):
            p.mkdir(parents=True, exist_ok=True)
