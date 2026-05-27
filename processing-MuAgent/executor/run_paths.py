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
            plan_summary.md             ← P2 output
            plan_review.md              ← `executor plan-review`
        post_run/                   materials produced DURING or AFTER execution
          summary/
            resolution_summary.md       ← S7 approval helper
            qc_summary.md               ← manifest rule
            run_manifest.json           ← manifest rule (handoff artifact)
            layout.json                 ← layout.finalize (manifest of deliverables)
          figures/                  user-facing figures only (QC + UMAP)
            s1_rna_qc_violin_{pre,post}.{png,pdf}
            s8_umap_{rna,atac}_by_leiden.{png,pdf}
          processed/
            processed.h5mu              (paired branch)
            rna_processed.h5ad          (separate / rna_only branches)
            atac_processed.h5ad         (separate / atac_only branches)
          notebooks/
            review_processed_h5mu.ipynb + .py
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
- Figures (QC + UMAP only) go to post_run/figures; no diagnostic figures anywhere.
- pre_run/ contains only files the user reviews before approving the plan.
- post_run/ contains in-run approval helpers (resolution_summary) + final results.
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

    # --- Deliverable phase roots ------------------------------------------
    @property
    def deliv_pre_run(self) -> Path:
        return self.deliverables / "pre_run"

    @property
    def deliv_post_run(self) -> Path:
        return self.deliverables / "post_run"

    # --- Deliverable sub-directories --------------------------------------
    @property
    def deliv_config(self) -> Path:
        """pre_run/config/ — user-supplied config + biological context."""
        return self.deliv_pre_run / "config"

    @property
    def deliv_pre_summary(self) -> Path:
        """pre_run/summary/ — P1 / P2 / plan_review materials for approval."""
        return self.deliv_pre_run / "summary"

    @property
    def deliv_post_summary(self) -> Path:
        """post_run/summary/ — in-run helpers + final run summaries."""
        return self.deliv_post_run / "summary"

    @property
    def deliv_figures(self) -> Path:
        return self.deliv_post_run / "figures"

    @property
    def deliv_processed(self) -> Path:
        return self.deliv_post_run / "processed"

    @property
    def deliv_notebooks(self) -> Path:
        return self.deliv_post_run / "notebooks"

    # --- Canonical user-facing files --------------------------------------
    # NOTE: these are deliverables (written directly), NOT internal copies.
    @property
    def run_yaml(self) -> Path:
        return self.deliv_config / "run.yaml"

    @property
    def biological_context_md(self) -> Path:
        return self.deliv_config / "biological_context.md"

    # pre_run/summary/*
    @property
    def context_summary_md(self) -> Path:
        return self.deliv_pre_summary / "context_summary.md"

    @property
    def plan_summary_md(self) -> Path:
        return self.deliv_pre_summary / "plan_summary.md"

    @property
    def plan_review_md(self) -> Path:
        return self.deliv_pre_summary / "plan_review.md"

    # post_run/summary/*
    @property
    def resolution_summary_md(self) -> Path:
        return self.deliv_post_summary / "resolution_summary.md"

    @property
    def qc_summary_pre_dimred_md(self) -> Path:
        """post_run/summary/qc_summary_pre_dimred.md — early QC review before S4/S5."""
        return self.deliv_post_summary / "qc_summary_pre_dimred.md"

    @property
    def qc_summary_md(self) -> Path:
        return self.deliv_post_summary / "qc_summary.md"

    @property
    def run_manifest_json(self) -> Path:
        return self.deliv_post_summary / "run_manifest.json"

    @property
    def layout_json(self) -> Path:
        return self.deliv_post_summary / "layout.json"

    # post_run/processed/*
    @property
    def processed_h5mu(self) -> Path:
        return self.deliv_processed / "processed.h5mu"

    @property
    def rna_processed_h5ad(self) -> Path:
        return self.deliv_processed / "rna_processed.h5ad"

    @property
    def atac_processed_h5ad(self) -> Path:
        return self.deliv_processed / "atac_processed.h5ad"

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

    def deliv_figure(self, stem: str, *, ext: str = "png") -> Path:
        return self.deliv_figures / f"{stem}.{ext}"

    # --- Bootstrap ---------------------------------------------------------
    def ensure(self) -> None:
        """Create both the internal and deliverables scaffold (idempotent)."""
        for p in (
            self.internal,
            self.artifacts, self.proposals, self.checkpoints,
            self.deliverables,
            self.deliv_pre_run, self.deliv_config, self.deliv_pre_summary,
            self.deliv_post_run, self.deliv_post_summary,
            self.deliv_figures, self.deliv_processed, self.deliv_notebooks,
        ):
            p.mkdir(parents=True, exist_ok=True)
