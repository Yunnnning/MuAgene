"""P2 preprocessing plan assembler.

Pulls approved P1 context + S0 validation and produces preprocessing_plan.json
with {value, source, rationale, confidence} per parameter.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import hashing as _h
from .hpc import load_execution_settings


# Per-branch stage set. Keys are stage IDs; values are the stages that should
# appear under `plan["stages"]` for that branch. See also: workflow_branch ∈
# {paired, separate, rna_only, atac_only}.
_STAGES_BY_BRANCH = {
    "paired":    {"s1a_ambient", "s1_rna_qc", "s2_atac_qc", "s3_doublets", "s4_rna_norm",
                   "s5_atac_lsi", "s6_dimred", "s7_clustering", "s8_umap"},
    "separate":  {"s1a_ambient", "s1_rna_qc", "s2_atac_qc", "s3_doublets", "s4_rna_norm",
                   "s5_atac_lsi", "s6_dimred", "s7_clustering", "s8_umap"},
    "rna_only":  {"s1a_ambient", "s1_rna_qc", "s3_doublets", "s4_rna_norm",
                   "s6_dimred", "s7_clustering", "s8_umap"},
    "atac_only": {"s2_atac_qc", "s3_doublets", "s5_atac_lsi",
                   "s6_dimred", "s7_clustering", "s8_umap"},
}


def _stages_for_branch(branch: str) -> set[str]:
    """Return the stage IDs that should appear in `plan['stages']` for `branch`."""
    try:
        return _STAGES_BY_BRANCH[branch]
    except KeyError as exc:
        raise ValueError(
            f"Unknown workflow_branch={branch!r}; expected one of {sorted(_STAGES_BY_BRANCH)}."
        ) from exc


_AMBIENT_METHODS = frozenset({"auto", "none", "decontx", "soupx"})


def _default_ambient_method(
    *,
    sample_type: str,
    ingest: dict[str, Any] | None,
    study_goal: str | None,
    user_method: str | None,
) -> tuple[str, str, str, str]:
    """Return (method, source, rationale, confidence) for s1a_ambient.method.

    Correction is dataset- and goal-driven (10x ambient-RNA guidance), not
    cells-vs-nuclei. Default is ``auto`` on RNA branches; confirm at plan review.
    """
    if user_method is not None:
        m = str(user_method).strip().lower()
        if m not in _AMBIENT_METHODS:
            raise ValueError(
                f"s1a_ambient_method must be one of {sorted(_AMBIENT_METHODS)} — got {user_method!r}."
            )
        return (
            m,
            "user",
            f"Explicit run.yaml override (s1a_ambient_method={m}).",
            "high",
        )

    has_raw = bool((ingest or {}).get("has_raw_matrix"))
    rna_status = (ingest or {}).get("rna_filtered_status")
    auto_hint = (
        "SoupX when raw + filtered matrices are present"
        if has_raw or rna_status == "raw"
        else "DecontX on filtered counts when no raw matrix"
    )
    goal = (study_goal or "clustering_inference").strip().lower()
    nuclei_note = (
        " Sample type is nuclei: ambient RNA can still be elevated with debris "
        "(10x guidance); disable only if you have inspected the data and "
        "contamination is low."
        if sample_type == "nuclei"
        else ""
    )
    if goal == "rare_populations":
        return (
            "auto",
            "recommended",
            f"study_goal=rare_populations: ambient correction recommended to "
            f"protect rare types from soup contamination. S0 dispatch: {auto_hint}."
            f"{nuclei_note} Set method=none at plan review if web summary / markers "
            "show low background.",
            "high",
        )
    return (
        "auto",
        "default",
        f"Default auto-dispatch from S0 inputs ({auto_hint}). Confirm at plan "
        f"review: set method=none if major cell types are clear and ambient "
        f"contamination is low (10x: not every dataset needs correction)."
        f"{nuclei_note}",
        "medium",
    )


def assemble_plan(
    run_dir: Path | str,
    *,
    workflow_branch: str,
    sample_type: str = "unknown",
    study_goal: str | None = None,
    ingest: dict[str, Any] | None = None,
    s1a_ambient_method: str | None = None,
) -> dict[str, Any]:
    run_dir = Path(run_dir)

    def p(value, source, rationale, confidence):
        return {"value": value, "source": source, "rationale": rationale, "confidence": confidence}

    ambient_method, amb_src, amb_rat, amb_conf = _default_ambient_method(
        sample_type=sample_type,
        ingest=ingest,
        study_goal=study_goal,
        user_method=s1a_ambient_method,
    )

    # Sample-type-aware ceilings for RNA QC
    if sample_type == "nuclei":
        pct_mt_ceil = 10.0
        pct_mt_rat = "Nuclei sample: cytoplasmic mRNA largely absent, mito content expected low."
    else:
        pct_mt_ceil = 20.0
        pct_mt_rat = "Whole-cell / unknown sample: standard ceiling."

    stages: dict[str, Any] = {
        "s1a_ambient": {
            "parameters": {
                "method": p(ambient_method, amb_src, amb_rat, amb_conf),
                "max_contamination": p(0.5, "recommended",
                                        "Cap per-cell rho/contamination at this fraction; "
                                        "prevents pathological over-correction on noisy cells.", "medium"),
            }
        },
        "s1_rna_qc": {
            "parameters": {
                "k_mad": p(5.0, "default", "Project convention for symmetric MAD on log1p counts.", "high"),
                "pct_mt_k": p(3.0, "default", "MAD multiplier for mito upper bound.", "high"),
                "pct_mt_ceiling": p(pct_mt_ceil, "inferred", pct_mt_rat, "medium"),
                "pct_mt_floor": p(5.0, "recommended", "Floor for pct_mt ceiling; avoids overly permissive cap on pristine samples.", "medium"),
                "pct_ribo_max": p(50.0, "recommended",
                                    "Soft ceiling on pct_counts_ribo (Rps/Rpl/Mrps/Mrpl). "
                                    "Stressed/dying cells often exceed this; tissues with very high "
                                    "ribo expression (e.g. plasma cells) may need a higher value.", "medium"),
                "min_cells_per_gene": p(3, "default", "scanpy convention.", "high"),
                "min_counts_floor": p(500, "recommended", "Guard against empty droplets dragging MAD down.", "medium"),
            }
        },
        "s2_atac_qc": {
            "parameters": {
                "tss_enrichment_min": p(1.5, "recommended",
                                        "Minimum TSS enrichment; cells at or below are removed.", "high"),
                "tss_enrichment_max": p(50.0, "recommended",
                                         "Maximum TSS enrichment; very high values often indicate artifacts.",
                                         "medium"),
                "n_fragments_k_mad": p(5.0, "default", "Symmetric MAD on log fragments per cell.", "high"),
                "n_fragments_floor": p(500, "recommended", "Minimum fragments for a real cell.", "medium"),
                "nucleosome_signal_max": p(2.0, "recommended",
                                            "Upper bound on nucleosome signal (mono/nucleosome-free fragment "
                                            "ratio). Cells at or above are removed.", "medium"),
            }
        },
        "s3_doublets": {
            "parameters": {
                "scrublet_expected_rate": p("auto", "recommended",
                                              "If 'auto', the rate scales as min(0.10, 0.0008 * n_cells) "
                                              "to track 10x's empirical doublet curve (~0.8% per 1000 cells). "
                                              "Override with a float to force a fixed rate.", "high"),
                "atac_doublet_threshold": p(0.5, "recommended",
                                              "SnapATAC2 scrublet doublet-probability threshold.", "medium"),
                "removal_policy_recommendation": p(
                    "independent" if workflow_branch == "separate" else "union",
                    "derived" if workflow_branch == "separate" else "recommended",
                    ("separate branch: each modality's doublets removed independently; "
                     "no cross-modal reconciliation." if workflow_branch == "separate" else
                     "Paired multiome: union of RNA and ATAC doublet calls (remove if either detector flags)."),
                    "high",
                ),
                "study_goal": p(study_goal or "clustering_inference", "user" if study_goal else "default",
                                 "From run.yaml or fallback.", "high" if study_goal else "medium"),
            }
        },
        "s4_rna_norm": {
            "parameters": {
                "target_sum": p(1e4, "default", "scanpy convention.", "high"),
                "hvg_flavor": p("seurat_v3", "default", "scanpy-native; operates on raw counts layer.", "high"),
                "hvg_n_top_genes": p(2000, "default", "Cap; actual count min(2000, 0.1 * n_genes_after_qc).", "high"),
            }
        },
        "s5_atac_lsi": {
            "parameters": {
                "n_components": p(50, "default", "Standard LSI dimensionality.", "high"),
                "drop_first": p(True, "default", "First LSI component is depth-correlated.", "high"),
                "max_top_peaks": p(50000, "recommended", "Cap on feature selection.", "medium"),
            }
        },
        "s6_dimred": {
            "parameters": {
                "rna_n_pcs": p("auto", "recommended",
                                "If 'auto', n_pcs is chosen by elbow detection on the cumulative "
                                "explained-variance curve (knee with `chord_distance`). Override "
                                "with an int to force a fixed value.", "medium"),
                "rna_n_pcs_max": p(50, "default",
                                     "Upper cap for the auto-elbow search.", "high"),
                "rna_scale": p(True, "recommended",
                                 "Apply sc.pp.scale(max_value=10) before PCA — the scanpy-standard "
                                 "preprocessing path. Disable to PCA the unscaled log-normalized data.", "high"),
                "n_neighbors": p(15, "default", "scanpy convention.", "high"),
            }
        },
        "s7_clustering": {
            "parameters": {
                # Backwards-compat single grid (used if per-modality grids are absent).
                "leiden_resolution_grid": p([0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2], "default",
                                             "Project standard grid (legacy; superseded by per-modality grids).", "high"),
                # Per-modality grids: ATAC is shifted to a lower range per the atac_tilt=lower policy.
                "leiden_resolution_grid_rna":  p([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2], "default",
                                                  "RNA grid; tolerates finer granularity (0.4 dropped; 1.1 added to close the 1.0→1.2 gap).", "high"),
                "leiden_resolution_grid_atac": p([0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], "default",
                                                  "ATAC grid; shifted lower to match atac_tilt=lower and avoid over-fragmentation.", "high"),
                "seeds": p([0, 1, 2], "default", "Three seeds for stability ARI.", "high"),
                "stability_floor": p(0.85, "default", "Minimum seed-pairwise ARI for stable region.", "medium"),
                "rna_tilt": p("higher", "default",
                               "RNA tolerates finer granularity per project policy.", "high"),
                "atac_tilt": p("lower", "default",
                               "ATAC prefers broader clusters to avoid over-fragmentation.", "high"),
            }
        },
        "s8_umap": {
            "parameters": {
                "min_dist": p(0.5, "default", "scanpy UMAP default.", "high"),
                "spread": p(1.0, "default", "scanpy UMAP default.", "high"),
                "random_state": p(42, "user", "Run seed from config.", "high"),
            }
        },
    }

    # Filter stages by branch — single-modality branches drop the irrelevant
    # per-modality stages (e.g. rna_only drops s2_atac_qc + s5_atac_lsi).
    keep = _stages_for_branch(workflow_branch)
    stages = {k: v for k, v in stages.items() if k in keep}

    plan: dict[str, Any] = {
        "workflow_branch": workflow_branch,
        "context_ref": "artifacts/p1_context/context_extraction.json",
        "ingest_ref": "artifacts/s0_ingest/validation_report.json",
        "execution": load_execution_settings(run_dir),
        "stages": stages,
        "assumptions": [
            f"sample_type={sample_type}",
            f"workflow_branch={workflow_branch}",
        ],
        "warnings": [],
    }
    return plan


def write_plan(run_dir: Path | str, plan: dict[str, Any]) -> tuple[Path, str]:
    from .run_paths import RunPaths
    out = RunPaths(Path(run_dir)).artifact("p2_plan", "preprocessing_plan.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(plan, indent=2, sort_keys=True, default=str)
    out.write_text(payload)
    phash = _h.sha256_bytes(payload.encode())
    return out, phash


def render_plan_appendix(plan: dict[str, Any]) -> str:
    """Per-stage parameter listing for the plan-review appendix."""
    lines: list[str] = [
        "## Appendix: full parameters",
        "",
        f"**Workflow branch:** `{plan['workflow_branch']}`",
        "",
    ]
    exec_block = plan.get("execution") or {}
    if exec_block:
        lines.extend(_render_execution_appendix(exec_block))
        lines.append("")
    for stage, body in plan["stages"].items():
        lines.append(f"### {stage}")
        for pname, pv in body["parameters"].items():
            lines.append(
                f"- **{pname}**: `{pv['value']}` — {pv['source']}/{pv['confidence']}"
            )
            lines.append(f"  - {pv['rationale']}")
        lines.append("")
    if plan.get("warnings"):
        lines.append("### Warnings")
        for w in plan["warnings"]:
            lines.append(f"- {w}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _render_execution_appendix(exec_block: dict[str, Any]) -> list[str]:
    mode = exec_block.get("mode", "local")
    settings = exec_block.get("settings") or {}
    lines = [
        "### Execution",
        f"- **mode**: `{mode}`",
    ]
    if exec_block.get("hpc_env_path"):
        lines.append(f"- **hpc_env**: `{exec_block['hpc_env_path']}`")
    if mode == "pbs":
        for key in ("pbs_queue", "pbs_project", "notify_email", "resources_scale", "conda_env"):
            val = settings.get(key)
            if val:
                lines.append(f"- **{key}**: `{val}`")
    elif mode == "slurm":
        for key in ("slurm_partition", "slurm_account", "notify_email", "resources_scale", "conda_env"):
            val = settings.get(key)
            if val:
                lines.append(f"- **{key}**: `{val}`")
    elif mode == "local":
        lines.append("- **note**: local mode; HPC vars not required unless S0 is retried on the cluster.")
    if exec_block.get("s0_policy"):
        lines.append(f"- **S0 policy**: {exec_block['s0_policy']}")
    return lines
