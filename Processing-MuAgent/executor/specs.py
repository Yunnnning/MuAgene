"""Per-stage job spec authoring — written by Processing-MuAgent, consumed by Execution-MuAgent.

Each spec is the formal statement of science intent for one pipeline stage:
what the biology needs (resources, inputs, outputs, science description) plus
the progress_timeout_hint that tells the monitor how long to wait before treating
the stage as hung. Platform mechanics (partition, account, container) come from
site.config and are never encoded here.

The spec I/O mirrors the Snakemake rule's DECLARED inputs/outputs — i.e. the durable
per-stage markers (validation_report.json, *_summary.json, calls.parquet, …), never the
deletable working h5ads. That is what survives cleanup and what the monitor should verify.
Branch membership and branch-aware inputs are derived from `pipeline` (the single source
of truth), so e.g. atac_only correctly gets s6_neighbors/s7_clustering specs whose inputs
key off the S5 spectral marker rather than the (absent) S4 RNA marker.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
from typing import Any

import yaml

from .pipeline import (
    PIPELINE_STAGE_ORDER,
    STAGES_BY_BRANCH,
    branch_has_atac,
    stages_for_branch,
)


_SCIENCE_DESCRIPTIONS: dict[str, str] = {
    "s0_ingest":        "Validate/ingest RNA+ATAC, assemble the preprocessing plan, and derive data-driven QC thresholds with per-metric removal previews and figures",
    "s1a_ambient":      "Correct ambient RNA contamination using SoupX or DecontX",
    "s1_rna_qc":        "Apply MAD-based QC thresholds to filter low-quality RNA cells",
    "s2_atac_qc":       "Apply TSS enrichment and nucleosome signal QC to filter low-quality ATAC cells",
    "s3_doublets":      "Detect and remove doublets using Scrublet (RNA) and SnapATAC2 (ATAC); fixed score thresholds",
    "s4_rna_norm":      "Normalize RNA counts and select highly variable genes",
    "s5_atac_spectral": "Compute ATAC spectral embedding from peak-by-cell matrix",
    "s6_neighbors":     "Build the neighbor graph: RNA PCA on RNA-bearing branches, ATAC KNN on the spectral embedding (X_spectral) when ATAC is present",
    "s7_clustering":    "Run Leiden clustering at fixed per-modality resolutions (RNA=0.7, ATAC=0.5)",
    "s8_umap":          "Project RNA and ATAC embeddings to 2D UMAP for visualization",
    "qc_handoff":        "Assemble the per-sample post-QC .h5mu + handoff manifest for Integration-MuAgent",
}

# Run-dir-relative input/output artifact paths per stage, mirroring each rule's DECLARED
# durable I/O (see module docstring). {run_dir}/{run} are templated, resolved at write time.
# Per-entry shape:
#   inputs           — always-present input markers
#   inputs_by_stage  — name -> (source_stage, path); included only when source_stage runs
#                      on this branch (branch-aware DAG edges, e.g. s3/s6). Ties directly to
#                      pipeline.STAGES_BY_BRANCH so spec inputs can never drift from membership.
#   outputs          — always-present declared outputs (the durable stage-done marker)
#   outputs_atac     — merged into outputs only on ATAC-bearing branches
_STAGE_IO: dict[str, dict[str, Any]] = {
    "s0_ingest": {
        "inputs": {
            "context": "{run_dir}/internal/artifacts/p1_context/context_extraction.json",
        },
        "outputs": {
            # validation_report.json is the durable S0 done-marker; rna_ingest.h5ad is a
            # deletable cache and is intentionally NOT a declared output.
            "validation_report": "{run_dir}/internal/artifacts/s0_ingest/validation_report.json",
            "plan":              "{run_dir}/internal/artifacts/p2_plan/preprocessing_plan.json",
            "qc_explore":        "{run_dir}/internal/artifacts/qc_explore/qc_explore.json",
        },
    },
    "s1a_ambient": {
        "inputs": {
            # Real DAG edge is S0's durable marker; rna_ingest.h5ad is a deletable cache
            # S1a reconstructs via io.load_rna_ingest when absent.
            "s0_marker": "{run_dir}/internal/artifacts/s0_ingest/validation_report.json",
        },
        "outputs": {
            # Sole declared output + durable stage-done marker (rna_decontaminated.h5ad is
            # a deletable working file read by S1 by path).
            "summary": "{run_dir}/internal/artifacts/s1a_ambient/summary.json",
        },
    },
    "s1_rna_qc": {
        "inputs": {
            # Durable S1a marker, not the deletable rna_decontaminated.h5ad.
            "s1a_marker": "{run_dir}/internal/artifacts/s1a_ambient/summary.json",
        },
        "outputs": {
            "qc_summary": "{run_dir}/internal/artifacts/s1_rna_qc/qc_summary.json",
        },
    },
    "s2_atac_qc": {
        "inputs": {
            "plan": "{run_dir}/internal/artifacts/p2_plan/preprocessing_plan.json",
        },
        "outputs": {
            "qc_summary": "{run_dir}/internal/artifacts/s2_atac_qc/qc_summary.json",
        },
    },
    "s3_doublets": {
        "inputs_by_stage": {
            # Durable upstream QC markers; present per branch (RNA branches get rna, ATAC
            # branches get atac), exactly mirroring s3_doublets.smk's _s3_inputs.
            "rna_qc_marker":  ("s1_rna_qc",  "{run_dir}/internal/artifacts/s1_rna_qc/qc_summary.json"),
            "atac_qc_marker": ("s2_atac_qc", "{run_dir}/internal/artifacts/s2_atac_qc/qc_summary.json"),
        },
        "outputs": {
            # Durable stage-done marker + DAG edge to qc_handoff. The post-doublet h5ads
            # are transient working files (deleted by qc_handoff), not declared.
            "calls": "{run_dir}/internal/artifacts/s3_doublets/calls.parquet",
        },
    },
    "s4_rna_norm": {
        "inputs": {
            "post_qc_h5mu": "{run_dir}/deliverables/qc/post_qc_{run}.h5mu",
        },
        "outputs": {
            # Durable marker, not the deletable rna_norm.h5ad (read by S6 by path).
            "summary": "{run_dir}/internal/artifacts/s4_rna_norm/norm_summary.json",
        },
    },
    "s5_atac_spectral": {
        "inputs": {
            "post_qc_h5mu": "{run_dir}/deliverables/qc/post_qc_{run}.h5mu",
        },
        "outputs": {
            "summary": "{run_dir}/internal/artifacts/s5_atac_spectral/spectral_summary.json",
        },
    },
    "s6_neighbors": {
        "inputs_by_stage": {
            # Branch-aware durable upstream markers, mirroring s6_neighbors.smk's _s6_inputs:
            # RNA branches depend on the S4 norm marker, ATAC branches on the S5 spectral
            # marker. On atac_only only the spectral marker is advertised (S4 does not run).
            "rna_norm_marker":      ("s4_rna_norm",      "{run_dir}/internal/artifacts/s4_rna_norm/norm_summary.json"),
            "atac_spectral_marker": ("s5_atac_spectral", "{run_dir}/internal/artifacts/s5_atac_spectral/spectral_summary.json"),
        },
        "outputs": {
            # Durable marker, not the deletable rna_neighbors.h5ad (read by S7 by path).
            "summary": "{run_dir}/internal/artifacts/s6_neighbors/neighbors_summary.json",
        },
    },
    "s7_clustering": {
        "inputs": {
            # Durable S6 marker, not the deletable rna_neighbors.h5ad.
            "s6_marker": "{run_dir}/internal/artifacts/s6_neighbors/neighbors_summary.json",
        },
        "outputs": {
            # Durable marker, not the deletable rna_clustered.h5ad / atac_leiden_labels.parquet.
            "summary": "{run_dir}/internal/artifacts/s7_clustering/clustering_summary.json",
        },
    },
    "s8_umap": {
        "inputs": {
            # Durable S7 marker, not the deletable rna_clustered.h5ad.
            "s7_marker": "{run_dir}/internal/artifacts/s7_clustering/clustering_summary.json",
        },
        "outputs": {
            "sentinel": "{run_dir}/internal/artifacts/s8_umap/s8_done.txt",
        },
    },
    "qc_handoff": {
        "inputs": {
            # Durable S3 marker (reads the transient post-doublet h5ads by path).
            "calls": "{run_dir}/internal/artifacts/s3_doublets/calls.parquet",
        },
        "outputs": {
            "post_qc_h5mu": "{run_dir}/deliverables/qc/post_qc_{run}.h5mu",
            "post_qc_manifest": "{run_dir}/deliverables/qc/post_qc_manifest.json",
        },
        # Merged into outputs on ATAC-bearing branches only (matches qc_handoff.smk).
        "outputs_atac": {
            "peaks_bed": "{run_dir}/deliverables/qc/peaks_{run}.bed",
        },
    },
}


def _spec_stages(branch: str) -> list[str]:
    """Stages that get a per-stage spec for `branch`, in workflow order.

    Derived from pipeline membership (the single source of truth) plus the two
    non-plan stages that always bracket a cluster run: the s0_ingest planning compute
    and the qc_handoff Integration bundle. Because this reads STAGES_BY_BRANCH directly,
    specs can never drift from the DAG (e.g. atac_only includes s6_neighbors/s7_clustering).
    """
    members = stages_for_branch(branch)
    ordered = [s for s in PIPELINE_STAGE_ORDER if s in members]
    return ["s0_ingest", *ordered, "qc_handoff"]


def _load_resources_smk() -> Any:
    """Load workflow/resources.smk as a Python module.

    resources.smk is pure Python with a non-standard extension. We use
    SourceFileLoader directly because spec_from_file_location may not
    resolve a loader for the .smk extension on all platforms.
    """
    import importlib.machinery
    import types
    smk_path = Path(__file__).resolve().parent.parent / "workflow" / "resources.smk"
    if not smk_path.exists():
        raise ImportError(f"resources.smk not found at {smk_path}")
    loader = importlib.machinery.SourceFileLoader("resources_smk", str(smk_path))
    mod = types.ModuleType("resources_smk")
    loader.exec_module(mod)
    return mod


def _resolve_io(io: dict[str, str], run_dir: str) -> dict[str, str]:
    run_name = Path(run_dir).name
    return {k: v.format(run_dir=run_dir, run=run_name) for k, v in io.items()}


def _stage_io_for_branch(stage: str, branch: str, run_dir: str) -> tuple[dict[str, str], dict[str, str]]:
    """Resolve a stage's (inputs, outputs) for `branch`, honouring branch-aware edges.

    inputs_by_stage edges are included only when their source stage runs on this branch
    (so they track pipeline membership); outputs_atac is merged on ATAC-bearing branches.
    """
    io = _STAGE_IO.get(stage, {})
    members = stages_for_branch(branch) if branch in STAGES_BY_BRANCH else set()

    inputs = dict(io.get("inputs", {}))
    for name, (src_stage, path) in io.get("inputs_by_stage", {}).items():
        if src_stage in members:
            inputs[name] = path

    outputs = dict(io.get("outputs", {}))
    if branch_has_atac(branch):
        outputs.update(io.get("outputs_atac", {}))

    return _resolve_io(inputs, run_dir), _resolve_io(outputs, run_dir)


def build_stage_spec(
    stage: str,
    run_dir: str | Path,
    *,
    resources: dict[str, int],
    runtime_min: int,
    progress_timeout_hint: int,
    branch: str,
) -> dict[str, Any]:
    """Build a spec dict for one stage."""
    run_dir_s = str(Path(run_dir).resolve())
    inputs, outputs = _stage_io_for_branch(stage, branch, run_dir_s)
    return {
        "schema_version": "1",
        "stage": stage,
        "science_description": _SCIENCE_DESCRIPTIONS.get(stage, stage),
        "resources": {
            "cpus": resources["cpus"],
            "mem_mb": resources["mem_mb"],
            "walltime_min": runtime_min,
        },
        "inputs":  inputs,
        "outputs": outputs,
        "progress_timeout_hint": progress_timeout_hint,
    }


def write_stage_specs(run_dir: Path | str, branch: str) -> list[Path]:
    """Write per-stage metadata YAMLs for all active stages in the given workflow branch.

    Resources and progress hints are loaded live from workflow/resources.smk so
    PMA_RESOURCES_SCALE is honoured at write time. Returns the list of written
    paths. Safe to call multiple times — existing files are overwritten.
    """
    from .run_paths import RunPaths

    res_mod = _load_resources_smk()
    resources: dict[str, dict[str, int]] = res_mod.RESOURCES
    runtime: dict[str, int] = res_mod.RUNTIME
    timeout_hint: dict[str, int] = res_mod.PROGRESS_TIMEOUT_HINT

    paths = RunPaths(Path(run_dir))
    paths.stage_meta_dir.mkdir(parents=True, exist_ok=True)

    stages = _spec_stages(branch)
    written: list[Path] = []
    for stage in stages:
        if stage not in resources:
            continue
        spec = build_stage_spec(
            stage,
            run_dir,
            resources=resources[stage],
            runtime_min=runtime[stage],
            progress_timeout_hint=timeout_hint.get(stage, 20),
            branch=branch,
        )
        out = paths.stage_meta(stage)
        out.write_text(yaml.safe_dump(spec, default_flow_style=False, sort_keys=False))
        written.append(out)
    # Prune orphan per-stage specs: any <stage>.yaml not part of the current branch
    # (left behind by a renamed stage — e.g. s_handoff -> qc_handoff — or a branch
    # change) would otherwise be read by Execution-MuAgent's per-stage output
    # verification, causing a stage to go unverified or be checked against stale
    # paths. head_job.yaml is written separately by write_head_job_spec and is always
    # preserved.
    keep = {p.name for p in written} | {"head_job.yaml"}
    for p in paths.stage_meta_dir.glob("*.yaml"):
        if p.name not in keep:
            p.unlink()
    return written


def write_head_job_spec(run_dir: Path | str, target: str) -> Path:
    """Write internal/stage_meta/head_job.yaml — the submission vehicle for the Snakemake head-job.

    Resources match the runner.slurm defaults (1 CPU, 4 GB, 24 h).
    The target field records the Snakemake target so execute-spec can pass it
    through to launch_runner.sh via PMA_TARGET.
    """
    from .run_paths import RunPaths

    run_dir_path = Path(run_dir).resolve()
    paths = RunPaths(run_dir_path)
    paths.stage_meta_dir.mkdir(parents=True, exist_ok=True)

    # Populate the head-job's declared outputs from the target stage so the monitor
    # can verify a clean exit (emit stage_output_verified) even during the planning
    # S0 submission, when per-stage specs do not exist yet. Only single-stage
    # `*_execute` targets resolve to a known output set; `*_propose` / `all`
    # (multi-stage phases) leave outputs empty — those phases already have per-stage
    # specs that cover verification, so the head-job stays a pure orchestrator there.
    from . import provenance
    branch = provenance.current_branch(str(paths.parameters_yaml))

    stage_stem = target.removesuffix("_execute")
    if stage_stem != target:
        _inputs, outputs = _stage_io_for_branch(stage_stem, branch, str(run_dir_path))
    else:
        outputs = {}

    # GPU preflight gating: tell Execution whether this run's stages actually include a
    # GPU-capable one (single source = _GPU_CAPABLE in resources.smk). execute-spec skips
    # the GPU env reconcile/pull when this is False, so a device=gpu run with no GPU
    # consumer (preprocessing today — _GPU_CAPABLE is empty) never pulls the container for
    # nothing. Per target: a single <stage>_execute → that stage; a multi-stage target
    # (all / *_propose) → the whole branch DAG (planning stages are CPU).
    res_mod = _load_resources_smk()
    gpu_capable: set[str] = getattr(res_mod, "_GPU_CAPABLE", set())
    if stage_stem != target:
        target_stages = {stage_stem}
    else:
        target_stages = set(_spec_stages(branch))
    gpu_stages_present = bool(gpu_capable & target_stages)

    spec: dict[str, Any] = {
        "schema_version": "1",
        "stage": "head_job",
        "science_description": "Snakemake orchestrator — submits and monitors all per-stage child jobs",
        "resources": {
            "cpus": 1,
            "mem_mb": 4000,
            "walltime_min": 1440,  # 24 h; matches runner.slurm defaults
        },
        "inputs": {
            "config": str(run_dir_path / "deliverables" / "plan" / "config" / "run.yaml"),
        },
        "outputs": outputs,
        "progress_timeout_hint": 120,  # 2 h silence on a Snakemake orchestrator is suspicious
        "snakemake_target": target,
        # True iff this run's stages include a GPU-capable one (_GPU_CAPABLE ∩ stages).
        # Execution's execute-spec gates the GPU env preflight on this — no GPU consumer,
        # no container pull. Empty _GPU_CAPABLE (preprocessing) → always False today.
        "gpu_stages_present": gpu_stages_present,
    }
    out = paths.stage_meta("head_job")
    out.write_text(yaml.safe_dump(spec, default_flow_style=False, sort_keys=False))
    return out
