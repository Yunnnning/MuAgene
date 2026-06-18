"""Per-stage job spec authoring — written by Processing-MuAgent, consumed by Execution-MuAgent.

Each spec is the formal statement of science intent for one pipeline stage:
what the biology needs (resources, inputs, outputs, science description) plus
the progress_timeout_hint that tells the monitor how long to wait before treating
the stage as hung. Platform mechanics (partition, account, container) come from
site.config and are never encoded here.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
from typing import Any

import yaml


_SCIENCE_DESCRIPTIONS: dict[str, str] = {
    "s0_ingest":        "Validate/ingest RNA+ATAC, assemble the preprocessing plan, and derive data-driven QC thresholds with per-metric removal previews and figures",
    "s1a_ambient":      "Correct ambient RNA contamination using SoupX or DecontX",
    "s1_rna_qc":        "Apply MAD-based QC thresholds to filter low-quality RNA cells",
    "s2_atac_qc":       "Apply TSS enrichment and nucleosome signal QC to filter low-quality ATAC cells",
    "s3_doublets":      "Detect and remove doublets using Scrublet (RNA) and SnapATAC2 (ATAC); fixed score thresholds",
    "s4_rna_norm":      "Normalize RNA counts and select highly variable genes",
    "s5_atac_spectral": "Compute ATAC spectral embedding from peak-by-cell matrix",
    "s6_neighbors":     "Build RNA PCA and shared nearest-neighbor graph",
    "s7_clustering":    "Run Leiden clustering at fixed per-modality resolutions (RNA=0.7, ATAC=0.5)",
    "s8_umap":          "Project RNA and ATAC embeddings to 2D UMAP for visualization",
    "s_handoff":        "Assemble the per-sample post-QC .h5mu + handoff manifest for Integration-MuAgent",
}

# Run-dir-relative input/output artifact paths per stage.
# Uses {run_dir} as a template token; resolved to absolute paths at write time.
_STAGE_IO: dict[str, dict[str, dict[str, str]]] = {
    "s0_ingest": {
        "inputs": {
            "context": "{run_dir}/internal/artifacts/p1_context/context_extraction.json",
        },
        "outputs": {
            "rna_h5ad":          "{run_dir}/internal/artifacts/s0_ingest/rna_ingest.h5ad",
            "validation_report": "{run_dir}/internal/artifacts/s0_ingest/validation_report.json",
            "plan":              "{run_dir}/internal/artifacts/p2_plan/preprocessing_plan.json",
            "qc_explore":        "{run_dir}/internal/artifacts/qc_explore/qc_explore.json",
        },
    },
    "s1a_ambient": {
        "inputs": {
            "rna_h5ad": "{run_dir}/internal/artifacts/s0_ingest/rna_ingest.h5ad",
        },
        "outputs": {
            "rna_h5ad": "{run_dir}/internal/artifacts/s1a_ambient/rna_decontaminated.h5ad",
        },
    },
    "s1_rna_qc": {
        "inputs": {
            "rna_h5ad": "{run_dir}/internal/artifacts/s1a_ambient/rna_decontaminated.h5ad",
        },
        "outputs": {
            "qc_summary_json": "{run_dir}/internal/artifacts/s1_rna_qc/qc_summary.json",
        },
    },
    "s2_atac_qc": {
        "inputs": {
            "plan": "{run_dir}/internal/artifacts/p2_plan/preprocessing_plan.json",
        },
        "outputs": {
            "qc_summary_json": "{run_dir}/internal/artifacts/s2_atac_qc/qc_summary.json",
        },
    },
    "s3_doublets": {
        "inputs": {
            "rna_h5ad":  "{run_dir}/internal/artifacts/s1_rna_qc/rna_qc.h5ad",
            "atac_h5ad": "{run_dir}/internal/artifacts/s2_atac_qc/atac_qc.h5ad",
        },
        "outputs": {
            "rna_post":  "{run_dir}/internal/artifacts/s3_doublets/rna_post_doublet.h5ad",
            "atac_post": "{run_dir}/internal/artifacts/s3_doublets/atac_post_doublet.h5ad",
        },
    },
    "s4_rna_norm": {
        "inputs": {
            "rna_h5ad": "{run_dir}/internal/artifacts/s3_doublets/rna_post_doublet.h5ad",
        },
        "outputs": {
            "rna_h5ad": "{run_dir}/internal/artifacts/s4_rna_norm/rna_norm.h5ad",
        },
    },
    "s5_atac_spectral": {
        "inputs": {
            "atac_h5ad": "{run_dir}/internal/artifacts/s3_doublets/atac_post_doublet.h5ad",
        },
        "outputs": {
            "summary": "{run_dir}/internal/artifacts/s5_atac_spectral/spectral_summary.json",
        },
    },
    "s6_neighbors": {
        "inputs": {
            "rna_h5ad": "{run_dir}/internal/artifacts/s4_rna_norm/rna_norm.h5ad",
        },
        "outputs": {
            "rna_neighbors": "{run_dir}/internal/artifacts/s6_neighbors/rna_neighbors.h5ad",
        },
    },
    "s7_clustering": {
        "inputs": {
            "rna_h5ad": "{run_dir}/internal/artifacts/s6_neighbors/rna_neighbors.h5ad",
        },
        "outputs": {
            "rna_clustered": "{run_dir}/internal/artifacts/s7_clustering/rna_clustered.h5ad",
        },
    },
    "s8_umap": {
        "inputs": {
            "rna_h5ad": "{run_dir}/internal/artifacts/s7_clustering/rna_clustered.h5ad",
        },
        "outputs": {
            "sentinel": "{run_dir}/internal/artifacts/s8_umap/s8_done.txt",
        },
    },
    "s_handoff": {
        "inputs": {
            "rna_post":  "{run_dir}/internal/artifacts/s3_doublets/rna_post_doublet.h5ad",
            "atac_post": "{run_dir}/internal/artifacts/s3_doublets/atac_post_doublet.h5ad",
        },
        "outputs": {
            "post_qc_manifest": "{run_dir}/deliverables/results/post_qc_manifest.json",
        },
    },
}

_BRANCH_STAGES: dict[str, list[str]] = {
    "paired": [
        "s0_ingest", "s1a_ambient", "s1_rna_qc", "s2_atac_qc", "s3_doublets",
        "s4_rna_norm", "s5_atac_spectral", "s6_neighbors", "s7_clustering", "s8_umap", "s_handoff",
    ],
    "separate": [
        "s0_ingest", "s1a_ambient", "s1_rna_qc", "s2_atac_qc", "s3_doublets",
        "s4_rna_norm", "s5_atac_spectral", "s6_neighbors", "s7_clustering", "s8_umap", "s_handoff",
    ],
    "rna_only": [
        "s0_ingest", "s1a_ambient", "s1_rna_qc", "s3_doublets",
        "s4_rna_norm", "s6_neighbors", "s7_clustering", "s8_umap", "s_handoff",
    ],
    "atac_only": [
        "s0_ingest", "s2_atac_qc", "s3_doublets",
        "s5_atac_spectral", "s8_umap", "s_handoff",
    ],
}


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
    return {k: v.format(run_dir=run_dir) for k, v in io.items()}


def build_stage_spec(
    stage: str,
    run_dir: str | Path,
    *,
    resources: dict[str, int],
    runtime_min: int,
    progress_timeout_hint: int,
) -> dict[str, Any]:
    """Build a spec dict for one stage."""
    run_dir_s = str(Path(run_dir).resolve())
    io = _STAGE_IO.get(stage, {"inputs": {}, "outputs": {}})
    return {
        "schema_version": "1",
        "stage": stage,
        "science_description": _SCIENCE_DESCRIPTIONS.get(stage, stage),
        "resources": {
            "cpus": resources["cpus"],
            "mem_mb": resources["mem_mb"],
            "walltime_min": runtime_min,
        },
        "inputs":  _resolve_io(io["inputs"],  run_dir_s),
        "outputs": _resolve_io(io["outputs"], run_dir_s),
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

    stages = _BRANCH_STAGES.get(branch, _BRANCH_STAGES["paired"])
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
        )
        out = paths.stage_meta(stage)
        out.write_text(yaml.safe_dump(spec, default_flow_style=False, sort_keys=False))
        written.append(out)
    return written


def write_head_job_spec(run_dir: Path | str, target: str) -> Path:
    """Write internal/stage_meta/head_job.yaml — the submission vehicle for the Snakemake head-job.

    Resources match the runner.slurm/runner.pbs defaults (1 CPU, 4 GB, 24 h).
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
    stage_stem = target.removesuffix("_execute")
    io = _STAGE_IO.get(stage_stem) if stage_stem != target else None
    outputs = _resolve_io(io["outputs"], str(run_dir_path)) if io else {}

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
        from . import provenance
        branch = provenance.current_branch(str(paths.parameters_yaml))
        target_stages = set(_BRANCH_STAGES.get(branch, _BRANCH_STAGES["paired"]))
    gpu_stages_present = bool(gpu_capable & target_stages)

    spec: dict[str, Any] = {
        "schema_version": "1",
        "stage": "head_job",
        "science_description": "Snakemake orchestrator — submits and monitors all per-stage child jobs",
        "resources": {
            "cpus": 1,
            "mem_mb": 4000,
            "walltime_min": 1440,  # 24 h; matches runner.slurm/runner.pbs defaults
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
