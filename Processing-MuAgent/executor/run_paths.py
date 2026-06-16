"""Single source of truth for canonical run-directory paths.

Top-level layout (direct-write; external inputs referenced via symlinks):

    <run_dir>/
      deliverables/                 user-facing; split by lifecycle phase
        plan/                       materials produced BEFORE preprocessing execution
          config/
            run.yaml                    ← canonical copy of the user config
            biological_context.md       ← canonical Biological Context Report
          summary/
            context_summary.md          ← P1 output
            plan_review_<run>.md        ← plan review gate (summary + parameter appendix)
            plan_summary_<run>.html     ← self-contained web version (figures embedded)
        figures/                    all pipeline figures (PNG + PDF), any stage
        checkpoints/                intermediate review reports (no figure files)
          qc_review/                  QC review checkpoint (summaries only)
            qc_review_<run>.md
            qc_summary_<run>.html
        results/                    final deliverables (data + manifest; no figures)
          review_processed_h5mu.{ipynb,py}
          processed.h5mu              (paired branch)
          rna_processed.h5ad          (separate / rna_only branches)
          atac_processed.h5ad         (separate / atac_only branches)
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
- External raw inputs (e.g. companion raw RNA matrix) are symlinked from S0, not copied.
- Derived stage artifacts have a single canonical path under internal/artifacts/.
- All figures live in deliverables/figures/; checkpoint dirs hold reports with embedded refs.
- results/ contains processed data, final notebook, and manifest only.
- plan/ contains only files the user reviews before approving the plan.
- figures/, checkpoints/, and results/ are created lazily when first written — not at init.
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

    def _existing_file(self, canonical: Path, legacy: Path) -> Path:
        """Prefer canonical path; fall back to legacy layout for older runs."""
        if canonical.is_file():
            return canonical
        if legacy.is_file():
            return legacy
        return canonical

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
        """Run-local Snakemake state directory."""
        return self.internal / "snakemake"

    # --- Deliverable phase roots ------------------------------------------
    @property
    def deliv_plan(self) -> Path:
        return self.deliverables / "plan"

    @property
    def deliv_figures(self) -> Path:
        """figures/ — canonical location for all pipeline figures (PNG + PDF)."""
        return self.deliverables / "figures"

    @property
    def deliv_checkpoints(self) -> Path:
        return self.deliverables / "checkpoints"

    @property
    def deliv_qc_review(self) -> Path:
        """checkpoints/qc_review/ — QC review reports (md/html only)."""
        return self.deliv_checkpoints / "qc_review"

    @property
    def deliv_results(self) -> Path:
        """results/ — flat final deliverables (data, notebook, manifest)."""
        return self.deliverables / "results"

    # --- plan sub-directories ---------------------------------------------
    @property
    def deliv_config(self) -> Path:
        """plan/config/ — user-supplied config + biological context."""
        return self.deliv_plan / "config"

    @property
    def deliv_plan_summary(self) -> Path:
        """plan/summary/ — P1 / P2 / plan_review materials for approval."""
        return self.deliv_plan / "summary"

    # --- Canonical user-facing files --------------------------------------
    @property
    def run_yaml(self) -> Path:
        return self._existing_file(
            self.deliv_config / "run.yaml",
            self.deliverables / "pre_run" / "config" / "run.yaml",
        )

    @property
    def biological_context_md(self) -> Path:
        canonical = self.deliv_config / "biological_context.md"
        return self._existing_file(
            canonical,
            self.deliverables / "pre_run" / "config" / "biological_context.md",
        )

    @property
    def hpc_env_sh(self) -> Path:
        """plan/config/hpc.env — source-able PMA_* exports for cluster runs."""
        canonical = self.deliv_config / "hpc.env"
        return self._existing_file(
            canonical,
            self.deliverables / "pre_run" / "config" / "hpc.env",
        )

    @property
    def site_config(self) -> Path:
        """plan/config/site.config — YAML platform description consumed by Execution-MuAgent."""
        canonical = self.deliv_config / "site.config"
        return self._existing_file(
            canonical,
            self.deliverables / "pre_run" / "config" / "site.config",
        )

    @property
    def stage_meta_dir(self) -> Path:
        return self.internal / "stage_meta"

    def stage_meta(self, stage: str) -> Path:
        return self.stage_meta_dir / f"{stage}.yaml"

    @property
    def context_summary_md(self) -> Path:
        return self.deliv_plan_summary / "context_summary.md"

    @property
    def plan_review_md(self) -> Path:
        return self.deliv_plan_summary / f"plan_review_{self.run_dir.name}.md"

    @property
    def plan_summary_html(self) -> Path:
        """Self-contained HTML plan review (figures embedded as data URIs)."""
        return self.deliv_plan_summary / f"plan_summary_{self.run_dir.name}.html"

    @property
    def preprocessing_plan(self) -> Path:
        """The assembled preprocessing plan JSON — single source of truth.

        Lives under the ``p2_plan`` artifact namespace for historical reasons:
        the standalone ``p2_plan`` rule was merged into ``s0_ingest`` (which now
        assembles the plan in-process), but the plan's artifact location was kept
        stable so existing runs and resume sessions keep working. Always go
        through this property rather than re-deriving ``artifact("p2_plan", ...)``
        so the location has exactly one definition.
        """
        return self.artifact("p2_plan", "preprocessing_plan.json")

    @property
    def plan_intro(self) -> Path:
        """Persisted plan-review intro paragraph (agent-authored prose).

        Stored so that every render path — the ``executor plan-review --intro``
        call, the ``plan_review_propose`` Snakemake rule, and any later
        re-render — reproduces the same intro in the run-scoped plan review
        markdown/HTML. Without this, a propose re-render would silently
        drop the intro because it is otherwise only a transient CLI argument.
        Co-located with the plan it introduces.
        """
        return self.artifact("p2_plan", "plan_intro.md")

    @property
    def validation_report(self) -> Path:
        """S0 ingest validation report JSON (pairing, counts, genome check)."""
        return self.artifact("s0_ingest", "validation_report.json")

    @property
    def qc_review_summary_md(self) -> Path:
        return self.deliv_qc_review / f"qc_review_{self.run_dir.name}.md"

    @property
    def qc_summary_html(self) -> Path:
        return self.deliv_qc_review / f"qc_summary_{self.run_dir.name}.html"

    @property
    def qc_summary_pre_dimred_md(self) -> Path:
        return self.qc_review_summary_md

    @property
    def run_manifest_json(self) -> Path:
        return self.deliv_results / "run_manifest.json"

    @property
    def layout_json(self) -> Path:
        return self.deliv_results / "layout.json"

    @property
    def processed_h5mu(self) -> Path:
        return self.deliv_results / "processed.h5mu"

    @property
    def rna_processed_h5ad(self) -> Path:
        return self.deliv_results / "rna_processed.h5ad"

    @property
    def atac_processed_h5ad(self) -> Path:
        return self.deliv_results / "atac_processed.h5ad"

    @property
    def review_notebook_ipynb(self) -> Path:
        return self.deliv_results / "review_processed_h5mu.ipynb"

    @property
    def review_notebook_py(self) -> Path:
        return self.deliv_results / "review_processed_h5mu.py"

    # --- Stage helpers -----------------------------------------------------
    def stage_dir(self, stage: str) -> Path:
        return self.artifacts / stage

    def artifact(self, stage: str, *parts: str) -> Path:
        return self.stage_dir(stage).joinpath(*parts)

    def proposal(self, stage: str, suffix: str = ".yaml") -> Path:
        return self.proposals / f"{stage}{suffix}"

    def awaiting_sentinel(self, stage: str) -> Path:
        return self.proposals / f"{stage}.awaiting_approval"

    def approved_sentinel(self, stage: str) -> Path:
        return self.checkpoints / f"{stage}.approved"

    def deliv_figures_path(self, stem: str, *, ext: str = "png") -> Path:
        return self.deliv_figures / f"{stem}.{ext}"

    def _legacy_figure_locations(self, stem: str, *, ext: str = "png") -> tuple[Path, ...]:
        d = self.deliverables
        return (
            d / "figure" / f"{stem}.{ext}",
            d / "checkpoint" / "qc_review" / "figures" / f"{stem}.{ext}",
            d / "checkpoint" / "qc_review" / f"{stem}.{ext}",
            d / "checkpoints" / "qc_review" / "figures" / f"{stem}.{ext}",
            d / "checkpoints" / "qc_review" / f"{stem}.{ext}",
            self.deliv_qc_review / "figures" / f"{stem}.{ext}",
            self.deliv_qc_review / f"{stem}.{ext}",
            d / "checkpoint" / "resolution_review" / f"{stem}.{ext}",
            self.deliv_checkpoints / "resolution_review" / f"{stem}.{ext}",
            d / "post_run" / f"{stem}.{ext}",
            self.deliv_results / f"{stem}.{ext}",
        )

    def resolve_figure(self, stem: str, *, ext: str = "png") -> Path:
        """Return figure path, falling back to legacy layouts for older runs."""
        canonical = self.deliv_figures_path(stem, ext=ext)
        if canonical.is_file():
            return canonical
        for legacy in self._legacy_figure_locations(stem, ext=ext):
            if legacy.is_file():
                return legacy
        return canonical

    def migrate_legacy_figures_to_central(self) -> list[Path]:
        """Move png/pdf from legacy deliverable subdirs into deliverables/figures/. Idempotent."""
        self.deliv_figures.mkdir(parents=True, exist_ok=True)
        moved: list[Path] = []
        src_dirs = (
            self.deliverables / "figure",
            self.deliverables / "checkpoint" / "qc_review" / "figures",
            self.deliverables / "checkpoint" / "qc_review",
            self.deliverables / "checkpoints" / "qc_review" / "figures",
            self.deliverables / "checkpoints" / "qc_review",
            self.deliverables / "checkpoint" / "resolution_review",
            self.deliv_checkpoints / "resolution_review",
            self.deliverables / "post_run",
            self.deliv_results,
        )
        for src_dir in src_dirs:
            if not src_dir.is_dir():
                continue
            for path in sorted(src_dir.iterdir()):
                if not path.is_file() or path.suffix not in {".png", ".pdf"}:
                    continue
                dest = self.deliv_figures / path.name
                if dest.exists():
                    path.unlink()
                else:
                    path.rename(dest)
                moved.append(dest)
        return moved

    # --- Bootstrap ---------------------------------------------------------
    def ensure(self) -> None:
        """Create internal + plan scaffold only (idempotent).

        figures/, checkpoints/, and results/ are created lazily when first written.
        """
        for p in (
            self.internal,
            self.artifacts, self.proposals, self.checkpoints, self.snakemake_workdir,
            self.deliverables,
            self.deliv_plan, self.deliv_config, self.deliv_plan_summary,
        ):
            p.mkdir(parents=True, exist_ok=True)
