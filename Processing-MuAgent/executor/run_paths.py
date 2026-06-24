"""Single source of truth for canonical run-directory paths.

Top-level layout (direct-write; external inputs referenced via symlinks):

    <run_dir>/
      deliverables/                 user-facing; split by lifecycle phase
        plan/                       materials produced BEFORE preprocessing execution
          config/
            run.yaml                    ← canonical copy of the user config
            biological_context.md       ← canonical Biological Context Report
          context_summary.md          ← P1 output
          plan_review_<run>.md        ← plan review gate (summary + parameter appendix)
          plan_summary_<run>.html     ← self-contained web version (figures embedded)
        figures/                    all pipeline figures (PNG + PDF), any stage
        qc/                         QC checkpoint: reports + post-QC Integration handoff
          qc_review_<run>.md
          qc_summary_<run>.html
          post_qc_<run>.h5mu         (after QC approval; all branches)
          peaks_<run>.bed            (after QC approval; ATAC branches with a peak set)
          post_qc_manifest.json
        results/                    final deliverables (data + manifest; no figures)
          review_processed_<run>.{ipynb,py}
          processed_<run>.h5mu        (paired branch)
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
- results/ contains S8 processed data, final notebook, and run manifest only.
- qc/ holds QC reports (before approval) and the post-QC Integration handoff (after approval).
- plan/ contains only files the user reviews before approving the plan.
- figures/, qc/, and results/ are created lazily when first written — not at init.
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
    def deliv_qc(self) -> Path:
        """qc/ — QC checkpoint: reports + post-QC Integration handoff (after approval)."""
        return self.deliverables / "qc"

    @property
    def deliv_results(self) -> Path:
        """results/ — flat final deliverables (data, notebook, manifest)."""
        return self.deliverables / "results"

    # --- plan sub-directories ---------------------------------------------
    @property
    def deliv_config(self) -> Path:
        """plan/config/ — user-supplied config + biological context."""
        return self.deliv_plan / "config"

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
        canonical = self.deliv_plan / "context_summary.md"
        return self._existing_file(canonical, self.deliv_plan / "summary" / "context_summary.md")

    @property
    def plan_review_md(self) -> Path:
        name = f"plan_review_{self.run_dir.name}.md"
        canonical = self.deliv_plan / name
        return self._existing_file(canonical, self.deliv_plan / "summary" / name)

    @property
    def plan_summary_html(self) -> Path:
        """Self-contained HTML plan review (figures embedded as data URIs)."""
        name = f"plan_summary_{self.run_dir.name}.html"
        canonical = self.deliv_plan / name
        return self._existing_file(canonical, self.deliv_plan / "summary" / name)

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
        return self.deliv_qc / f"qc_review_{self.run_dir.name}.md"

    @property
    def qc_summary_html(self) -> Path:
        return self.deliv_qc / f"qc_summary_{self.run_dir.name}.html"

    @property
    def post_qc_h5mu(self) -> Path:
        name = f"post_qc_{self.run_dir.name}.h5mu"
        return self._existing_file(self.deliv_qc / name, self.deliv_results / name)

    @property
    def post_qc_manifest_json(self) -> Path:
        return self._existing_file(
            self.deliv_qc / "post_qc_manifest.json",
            self.deliv_results / "post_qc_manifest.json",
        )

    @property
    def post_qc_peaks_bed(self) -> Path:
        """Per-sample peaks BED in the QC integration bundle (moved/copied by qc_handoff)."""
        return self.deliv_qc / f"peaks_{self.run_dir.name}.bed"

    def resolve_post_qc_peaks_bed(self) -> Path | None:
        """Canonical post-QC peaks BED for S5 / Integration (legacy fallbacks included)."""
        if self.post_qc_peaks_bed.is_file():
            return self.post_qc_peaks_bed
        man = self.post_qc_manifest_json
        if man.is_file():
            try:
                import json
                rel = (json.loads(man.read_text()).get("atac") or {}).get("peaks_bed")
                if rel:
                    p = self.run_dir / rel
                    if p.is_file():
                        return p
            except Exception:
                pass
        for name in ("peaks_macs3.bed", "peaks_arc.bed"):
            p = self.artifact("s2_atac_qc", name)
            if p.is_file():
                return p
        return None

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
        canonical = self.deliv_results / f"processed_{self.run_dir.name}.h5mu"
        return self._existing_file(canonical, self.deliv_results / "processed.h5mu")

    @property
    def rna_processed_h5ad(self) -> Path:
        return self.deliv_results / "rna_processed.h5ad"

    @property
    def atac_processed_h5ad(self) -> Path:
        return self.deliv_results / "atac_processed.h5ad"

    @property
    def review_notebook_ipynb(self) -> Path:
        canonical = self.deliv_results / f"review_processed_{self.run_dir.name}.ipynb"
        return self._existing_file(canonical, self.deliv_results / "review_processed_h5mu.ipynb")

    @property
    def review_notebook_py(self) -> Path:
        canonical = self.deliv_results / f"review_processed_{self.run_dir.name}.py"
        return self._existing_file(canonical, self.deliv_results / "review_processed_h5mu.py")

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
            d / "checkpoint" / "resolution_review" / f"{stem}.{ext}",
            d / "checkpoints" / "resolution_review" / f"{stem}.{ext}",
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
            self.deliverables / "checkpoint" / "resolution_review",
            self.deliverables / "checkpoints" / "resolution_review",
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

        figures/, qc/, and results/ are created lazily when first written.
        """
        for p in (
            self.internal,
            self.artifacts, self.proposals, self.checkpoints, self.snakemake_workdir,
            self.deliverables,
            self.deliv_plan, self.deliv_config,
        ):
            p.mkdir(parents=True, exist_ok=True)
