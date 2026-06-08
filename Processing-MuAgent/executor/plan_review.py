"""User-facing plan review: concise summary built from P1 context + S0 ingest + P2 plan."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import provenance


def _load_json(p: Path) -> dict:
    return json.loads(p.read_text()) if p.exists() else {}


_CONTEXT_FIELD_LABELS: dict[str, str] = {
    "organism": "organism",
    "tissue": "tissue",
    "assay_type": "assay",
    "genome_build": "reference genome",
    "dois": "publication DOIs",
    "modality_type": "dataset type",
    "sample_type": "sample type",
}


def _pairing_review_reason(pairing: dict[str, Any]) -> str:
    method = pairing.get("method", "")
    overlap = float(pairing.get("overlap", 0) or 0)
    labels = {
        "pairing.single_file_multiome": (
            "RNA and ATAC from the same Cell Ranger ARC output; "
            "barcodes linked via the combined matrix file."
        ),
        "pairing.exact_barcode_match": "RNA and ATAC barcodes match directly.",
        "pairing.prefix_suffix_normalized": (
            "RNA and ATAC barcodes matched after normalising prefix/suffix differences."
        ),
        "pairing.translation_table": "RNA and ATAC barcodes matched via a user-supplied translation table.",
        "pairing.rna_subset_of_atac": "RNA barcodes are largely a subset of ATAC barcodes.",
        "pairing.atac_subset_of_rna": "ATAC barcodes are largely a subset of RNA barcodes.",
        "pairing.rna_only_input": "Only RNA input was provided.",
        "pairing.atac_only_input": "Only ATAC input was provided.",
        "pairing.ambiguous_overlap": f"Barcode overlap is ambiguous ({overlap:.1%}); confirm pairing before continuing.",
        "pairing.no_match": "No reliable RNA–ATAC barcode overlap detected.",
    }
    if method in labels:
        return labels[method]
    if method:
        return f"Pairing established during ingest ({overlap:.1%} barcode overlap)."
    return "Pairing assessed during ingest."


def build_intro_context(run_dir: Path | str) -> dict[str, Any]:
    """Return the flat data dict the agent needs to write the intro paragraph."""
    from .run_paths import RunPaths
    paths = RunPaths(Path(run_dir))
    ctx = _load_json(paths.artifact("p1_context", "context_extraction.json"))
    ingest = _load_json(paths.artifact("s0_ingest", "validation_report.json"))
    plan = _load_json(paths.artifact("p2_plan", "preprocessing_plan.json"))

    def field(name: str) -> dict[str, Any]:
        return ctx.get("fields", {}).get(name, {})

    pairing = ingest.get("pairing", {})
    overlap_raw = pairing.get("overlap")
    study_goal = (
        plan.get("stages", {})
        .get("s3_doublets", {})
        .get("parameters", {})
        .get("study_goal", {})
        .get("value", "")
    )

    return {
        "organism":             field("organism").get("value", ""),
        "tissue":               field("tissue").get("value", ""),
        "assay_type":           field("assay_type").get("value", ""),
        "sample_type":          field("sample_type").get("value", ""),
        "genome_build":         field("genome_build").get("value", ""),
        "workflow_branch":      plan.get("workflow_branch", ingest.get("workflow_branch", "")),
        "study_goal":           study_goal,
        "rna_n_cells":          ingest.get("rna_n_cells"),
        "atac_n_barcodes":      ingest.get("atac_n_unique_barcodes"),
        "rna_raw_n_barcodes":   ingest.get("rna_raw_n_barcodes"),
        "has_raw_matrix":       bool(ingest.get("has_raw_matrix", False)),
        "pairing_status":       pairing.get("status", ""),
        "pairing_overlap":      float(overlap_raw) if overlap_raw is not None else None,
        "pairing_method":       pairing.get("method", ""),
        "pairing_confidence":   pairing.get("confidence", ""),
        "single_file_multiome": bool(ingest.get("single_file_multiome", False)),
        "rna_filtered_status":  ingest.get("rna_filtered_status", ""),
        "atac_barcodes_source": ingest.get("atac_barcodes_source", ""),
        "atac_fragments_file_n": ingest.get("atac_fragments_file_barcodes_n"),
        "pairing_ladder":       pairing.get("ladder", []),
    }


def build_summary(run_dir: Path | str) -> list[dict[str, Any]]:
    """Produce an ordered list of review items: {label, value, reason, certainty}."""
    from .run_paths import RunPaths
    paths = RunPaths(Path(run_dir))
    run_dir = paths.run_dir
    ctx = _load_json(paths.artifact("p1_context", "context_extraction.json"))
    ingest = _load_json(paths.artifact("s0_ingest", "validation_report.json"))
    plan = _load_json(paths.artifact("p2_plan", "preprocessing_plan.json"))

    def field(name: str) -> dict[str, Any]:
        return ctx.get("fields", {}).get(name, {})

    def param(stage: str, key: str) -> dict[str, Any]:
        return plan.get("stages", {}).get(stage, {}).get("parameters", {}).get(key, {})

    items: list[dict[str, Any]] = []

    # 1. Execution mode + HPC settings (from configure-execution / hpc.env)
    items.extend(_execution_review_items(run_dir, plan))

    # 2. Dataset type / modality
    modality = field("modality_type")
    items.append({
        "label": "Detected dataset type",
        "value": modality.get("value", "unknown"),
        "reason": modality.get("rationale", ""),
        "certainty": "certain" if modality.get("status") == "explicit" else "needs confirmation",
    })

    organism = field("organism")
    items.append({
        "label": "Organism",
        "value": organism.get("value", "unknown"),
        "reason": organism.get("rationale", ""),
        "certainty": "certain" if organism.get("status") == "explicit" else "needs confirmation",
    })

    genome = field("genome_build")
    items.append({
        "label": "Genome reference",
        "value": genome.get("value", "unknown"),
        "reason": genome.get("rationale", ""),
        "certainty": "certain" if genome.get("status") == "explicit" else "needs confirmation",
    })

    # 3. Pairing
    pairing = ingest.get("pairing", {})
    items.append({
        "label": "Pairing (RNA ↔ ATAC)",
        "value": f"{pairing.get('status', 'unknown')} (overlap={pairing.get('overlap', 0):.2%})",
        "reason": _pairing_review_reason(pairing),
        "certainty": "certain" if pairing.get("confidence") == "high" else "needs confirmation",
    })

    # 4. Ambient RNA correction (S1a) — confirm at plan review
    ambient_method = param("s1a_ambient", "method")
    rna_filtered_status = ingest.get("rna_filtered_status")
    has_raw = ingest.get("has_raw_matrix", False)
    if ambient_method.get("value") is not None:
        method_label = ambient_method.get("value", "auto")
        if method_label == "none":
            ambient_value = "Off — no ambient RNA correction"
            ambient_reason = (
                "Correction is disabled for this run. Turn it back on at plan review "
                "if background noise or contamination is a concern."
            )
        else:
            if method_label == "soupx":
                dispatch = "SoupX"
            elif method_label == "decontx":
                dispatch = "DecontX"
            elif has_raw and rna_filtered_status == "filtered":
                dispatch = "SoupX (filtered + raw RNA → auto)"
            else:
                dispatch = "DecontX (filtered RNA only → auto)"
            ambient_value = f"On — {dispatch}"
            ambient_reason = (
                "Available RNA inputs only pick the method: filtered alone uses "
                "DecontX; filtered plus raw uses SoupX. Whether to run ambient RNA "
                "correction depends on background noise and contamination — approve "
                "as-is if background RNA is a concern, or skip at plan review if "
                "contamination looks low after inspecting the data."
            )
        items.append({
            "label": "Ambient RNA correction (S1a)",
            "value": ambient_value,
            "reason": ambient_reason,
            "certainty": (
                "certain" if ambient_method.get("source") == "user" else "needs confirmation"
            ),
        })

    # 5. Key QC strategy
    sample_type = field("sample_type").get("value", "unknown")
    pct_mt_ceil = param("s1_rna_qc", "pct_mt_ceiling")
    pct_ribo_max = param("s1_rna_qc", "pct_ribo_max")
    tss_min = param("s2_atac_qc", "tss_enrichment_min")
    tss_max = param("s2_atac_qc", "tss_enrichment_max")
    nuc_max = param("s2_atac_qc", "nucleosome_signal_max")
    frip_min = param("s2_atac_qc", "frip_min")
    frip_val = frip_min.get("value", 0.0)
    # Peak source knowable at plan time: user-supplied path (highest priority) or ARC h5.
    # MACS3 is a runtime fallback — flag as conditional if neither is confirmed.
    has_user_peaks = bool(ingest.get("atac_peaks_path"))
    has_arc_peaks = bool(ingest.get("single_file_multiome"))
    if frip_val and float(frip_val) > 0:
        if has_user_peaks or has_arc_peaks:
            frip_note = f", FRiP >= {frip_val}"
        else:
            frip_note = f", FRiP >= {frip_val} (if peaks available at runtime)"
    else:
        frip_note = ""
    items.append({
        "label": "QC strategy",
        "value": (
            f"RNA: MAD on total_counts/n_genes, pct_mt ceiling={pct_mt_ceil.get('value', '?')}, "
            f"pct_ribo ceiling={pct_ribo_max.get('value', '?')} | "
            f"ATAC: TSS in ({tss_min.get('value', '?')}, {tss_max.get('value', '?')}), "
            f"MAD on log(n_fragments), nucleosome_signal<{nuc_max.get('value', '?')}{frip_note}"
        ),
        "reason": f"Sample type = {sample_type}; MAD thresholds adapt to the observed distribution.",
        "certainty": "certain",
    })

    # 5. Doublet policy — reconciliation is only meaningful when both detectors run.
    # Single-modality branches have one detector; the S3 item is informational, not
    # a checkpoint gate.
    branch = plan.get("workflow_branch", "paired")
    policy = param("s3_doublets", "removal_policy_recommendation")
    goal = param("s3_doublets", "study_goal")

    rna_doub_thr = param("s3_doublets", "rna_doublet_score_threshold")
    atac_doub_thr = (
        param("s3_doublets", "atac_doublet_probability_threshold")
        or param("s3_doublets", "atac_doublet_threshold")
        or param("s3_doublets", "atac_doublet_score_threshold")
    )

    def _paired_doublet_policy_detail(recommended: str) -> list[str]:
        lines = [
            "- paired-multiome policy: Each modality runs its own doublet detector. "
            "Cells flagged by **either** detector are removed (union). Detectors are "
            "prone to false negatives, so union minimises doublet contamination.",
            f"- applied policy: `{recommended}`",
        ]
        if rna_doub_thr.get("value") is not None:
            lines.append(
                f"- RNA Scrublet score threshold: `{rna_doub_thr.get('value')}` "
                "(cells with score above are flagged)"
            )
        if atac_doub_thr.get("value") is not None:
            lines.append(
                f"- ATAC SnapATAC2 probability threshold: `{atac_doub_thr.get('value')}` "
                "(cells with doublet_probability above are flagged)"
            )
        return lines

    if branch == "paired":
        policy_val = policy.get("value", "?")
        items.append({
            "label": "Doublet removal policy",
            "value": policy_val,
            "detail": _paired_doublet_policy_detail(policy_val),
            "reason": "Paired multiome always uses union: remove if either detector flags.",
            "certainty": "certain",
        })
    elif branch == "separate":
        detail = [
            "- separate branch: each modality's doublets removed independently "
            "(Scrublet for RNA, SnapATAC2 for ATAC).",
        ]
        if rna_doub_thr.get("value") is not None:
            detail.append(f"- RNA Scrublet score threshold: `{rna_doub_thr.get('value')}`")
        if atac_doub_thr.get("value") is not None:
            detail.append(f"- ATAC SnapATAC2 score threshold: `{atac_doub_thr.get('value')}`")
        items.append({
            "label": "Doublet removal policy",
            "value": "independent (per-modality; fixed score thresholds)",
            "detail": detail,
            "reason": ("separate branch: modalities are independent samples with disjoint barcodes; "
                       "each modality's doublets are removed by its own detector."),
            "certainty": "certain",
        })
    elif branch == "rna_only":
        rna_thr_note = ""
        if rna_doub_thr.get("value") is not None:
            rna_thr_note = f"; score threshold={rna_doub_thr.get('value')}"
        items.append({
            "label": "Doublet removal policy",
            "value": f"Scrublet only (fixed score threshold{rna_thr_note})",
            "reason": f"study_goal={goal.get('value', '?')}; single-detector branch — no reconciliation to confirm.",
            "certainty": "certain",
        })
    elif branch == "atac_only":
        atac_thr_note = ""
        if atac_doub_thr.get("value") is not None:
            atac_thr_note = f"; score threshold={atac_doub_thr.get('value')}"
        items.append({
            "label": "Doublet removal policy",
            "value": f"SnapATAC2 scrublet only (fixed score threshold{atac_thr_note})",
            "reason": f"study_goal={goal.get('value', '?')}; single-detector branch — no reconciliation to confirm.",
            "certainty": "certain",
        })

    # 6. Clustering
    grid = param("s7_clustering", "leiden_resolution_grid")
    rna_tilt = param("s7_clustering", "rna_tilt")
    atac_tilt = param("s7_clustering", "atac_tilt")
    items.append({
        "label": "Clustering strategy",
        "value": f"Leiden sweep {grid.get('value', [])} → stable-region knee (RNA={rna_tilt.get('value', '?')}, ATAC={atac_tilt.get('value', '?')})",
        "reason": "Per-modality resolution picked from a stability plateau, not a single-metric optimum.",
        "certainty": "certain",
    })

    # 7. Output location
    branch = plan.get("workflow_branch", "unknown")
    _s8 = paths.stage_dir("s8_umap")
    outputs = {
        "paired": f"{_s8}/processed.h5mu",
        "separate": f"{_s8}/rna_processed.h5ad + atac_processed.h5ad",
    }.get(branch, "(branch unknown)")
    branch_labels = {
        "paired": "Paired multiome run",
        "separate": "Separate RNA and ATAC samples",
        "rna_only": "RNA-only run",
        "atac_only": "ATAC-only run",
    }
    items.append({
        "label": "Output location",
        "value": outputs,
        "reason": (
            f"{branch_labels.get(branch, 'Run')}; "
            "a deliverable manifest is written alongside the processed object."
        ),
        "certainty": "certain",
    })

    # 8. Missing / uncertain
    missing = [k for k, v in ctx.get("fields", {}).items() if isinstance(v, dict) and v.get("status") == "missing"]
    low_conf = [k for k, v in ctx.get("fields", {}).items() if isinstance(v, dict) and v.get("confidence") == "low"]
    conflicts = ctx.get("conflicts", [])
    parts = []
    if missing:
        labels = [_CONTEXT_FIELD_LABELS.get(k, k.replace("_", " ")) for k in missing]
        parts.append(f"not provided: {', '.join(labels)}")
    if low_conf:
        labels = [_CONTEXT_FIELD_LABELS.get(k, k.replace("_", " ")) for k in low_conf]
        parts.append(f"low confidence: {', '.join(labels)}")
    if conflicts:
        parts.append(f"{len(conflicts)} unresolved conflict(s)")
    items.append({
        "label": "Missing / uncertain info",
        "value": "; ".join(parts) if parts else "none",
        "reason": "Gaps or uncertainty flagged during biological context intake.",
        "certainty": "certain" if not parts else "needs confirmation",
    })

    return items


def _execution_review_items(run_dir: Path, plan: dict[str, Any]) -> list[dict[str, Any]]:
    from .run_paths import RunPaths

    exec_block = plan.get("execution") or {}
    mode = exec_block.get("mode", "local")
    settings = exec_block.get("settings") or {}
    params = provenance.load(RunPaths(run_dir).parameters_yaml)
    mode_entry = params.get("execution.mode") or {}
    user_configured = isinstance(mode_entry, dict) and mode_entry.get("source") == "user"

    items: list[dict[str, Any]] = [{
        "label": "Execution mode",
        "value": mode,
        "reason": exec_block.get(
            "s0_policy",
            "How heavy stages (S0 onward) are dispatched: local machine vs PBS/SLURM.",
        ),
        "certainty": "certain" if user_configured else "needs confirmation",
    }]

    if mode == "local":
        hpc_path = exec_block.get("hpc_env_path")
        items.append({
            "label": "HPC configuration",
            "value": hpc_path or "none (local default)",
            "reason": (
                "Cluster vars not required for a fully local run. "
                "If S0 fails locally, the agent may configure HPC and retry ingest on the cluster."
            ),
            "certainty": "certain" if not hpc_path else "needs confirmation",
        })
        return items

    # PBS or SLURM — surface scheduler settings for review
    parts: list[str] = []
    missing: list[str] = []
    if mode == "pbs":
        for label, key in (("queue", "pbs_queue"), ("project", "pbs_project")):
            val = settings.get(key)
            if val:
                parts.append(f"{label}={val}")
            elif key == "pbs_queue":
                missing.append("pbs_queue")
        for label, key in (("scale", "resources_scale"), ("conda", "conda_env")):
            val = settings.get(key)
            if val:
                parts.append(f"{label}={val}")
    else:
        for label, key in (("partition", "slurm_partition"), ("account", "slurm_account")):
            val = settings.get(key)
            if val:
                parts.append(f"{label}={val}")
            elif key == "slurm_partition":
                missing.append("slurm_partition")
        for label, key in (("scale", "resources_scale"), ("conda", "conda_env")):
            val = settings.get(key)
            if val:
                parts.append(f"{label}={val}")

    items.append({
        "label": "HPC configuration",
        "value": ", ".join(parts) if parts else "not set",
        "reason": (
            "Cluster scheduler settings for this run. "
            + (f"Still needed: {', '.join(missing)}." if missing else "Recorded and ready for submit/resume.")
        ),
        "certainty": "needs confirmation" if missing else ("certain" if user_configured else "needs confirmation"),
    })
    return items


def _render_concise_section(items: list[dict[str, Any]]) -> str:
    lines: list[str] = []
    for it in items:
        tag = "✓" if it["certainty"] == "certain" else "?"
        lines.append(f"- **{it['label']}** [{tag} {it['certainty']}]")
        lines.append(f"  - value: `{it['value']}`")
        for detail_line in it.get("detail", []):
            lines.append(f"  {detail_line}")
        if it.get("reason"):
            lines.append(f"  - reason: {it['reason']}")
    return "\n".join(lines)


def render_merged_markdown(run_dir: Path | str, intro: str | None = None) -> str:
    """Full plan-review deliverable: concise summary + parameter appendix."""
    from .plan_assembler import render_plan_appendix
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    paths = RunPaths(run_dir)
    items = build_summary(run_dir)
    plan = _load_json(paths.artifact("p2_plan", "preprocessing_plan.json"))
    parts = ["# Preprocessing plan review", ""]
    if intro:
        parts += [intro.strip(), ""]
    parts += [
        "## Summary",
        "",
        _render_concise_section(items),
        "",
        "_Machine-readable plan: `internal/artifacts/p2_plan/preprocessing_plan.json` "
        "and `internal/parameters.yaml`. Full parameter listing in the appendix below._",
        "",
        "---",
        "",
    ]
    if plan:
        parts.append(render_plan_appendix(plan).rstrip())
        parts.append("")
    else:
        parts.append("_Appendix unavailable — `preprocessing_plan.json` not found._")
        parts.append("")
    return "\n".join(parts)


def write_summary(run_dir: Path | str, intro: str | None = None) -> Path:
    """Write merged plan-review markdown (summary + parameter appendix)."""
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    paths = RunPaths(run_dir)
    merged = render_merged_markdown(run_dir, intro=intro)
    paths.plan_review_md.parent.mkdir(parents=True, exist_ok=True)
    paths.plan_review_md.write_text(merged)
    return paths.plan_review_md
