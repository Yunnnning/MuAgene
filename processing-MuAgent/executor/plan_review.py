"""User-facing plan review: concise summary built from P1 context + S0 ingest + P2 plan."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import provenance


def _load_json(p: Path) -> dict:
    return json.loads(p.read_text()) if p.exists() else {}


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

    # 1. Dataset type / modality
    modality = field("modality_type")
    items.append({
        "label": "Detected dataset type",
        "value": modality.get("value", "unknown"),
        "reason": modality.get("rationale", ""),
        "certainty": "certain" if modality.get("status") == "explicit" else "needs confirmation",
    })

    # 2. Organism + genome build
    organism = field("organism")
    genome = field("genome_build")
    items.append({
        "label": "Organism / genome build",
        "value": f"{organism.get('value', 'unknown')} / {genome.get('value', 'unknown')}",
        "reason": (f"{organism.get('rationale', '')} | {genome.get('rationale', '')}").strip(" |"),
        "certainty": "certain" if organism.get("status") == "explicit" and genome.get("confidence") != "low" else "needs confirmation",
    })

    # 3. Pairing
    pairing = ingest.get("pairing", {})
    items.append({
        "label": "Pairing (RNA ↔ ATAC)",
        "value": f"{pairing.get('status', 'unknown')} (overlap={pairing.get('overlap', 0):.2%})",
        "reason": f"{pairing.get('method', '')} at confidence={pairing.get('confidence', '')}; shared barcodes={pairing.get('n_shared', 0)}.",
        "certainty": "certain" if pairing.get("confidence") == "high" else "needs confirmation",
    })

    # 4. Key QC strategy
    sample_type = field("sample_type").get("value", "unknown")
    pct_mt_ceil = param("s1_rna_qc", "pct_mt_ceiling")
    tss_min = param("s2_atac_qc", "tss_enrichment_min")
    items.append({
        "label": "QC strategy",
        "value": (
            f"RNA: MAD-based on total_counts/n_genes, pct_mt ceiling={pct_mt_ceil.get('value', '?')} | "
            f"ATAC: TSS≥{tss_min.get('value', '?')}, MAD on log(n_fragments)"
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
    if branch in ("paired", "separate"):
        items.append({
            "label": "Doublet removal policy",
            "value": f"{policy.get('value', '?')} (reconciling Scrublet RNA + ATAC detector; raw calls preserved)",
            "reason": f"study_goal={goal.get('value', '?')}; four-way overlap recorded; user confirms reconciliation at S3.",
            "certainty": "needs confirmation",
        })
    elif branch == "rna_only":
        items.append({
            "label": "Doublet removal policy",
            "value": "auto-applied (Scrublet only; raw calls preserved)",
            "reason": f"study_goal={goal.get('value', '?')}; single-detector branch — no reconciliation to confirm.",
            "certainty": "certain",
        })
    elif branch == "atac_only":
        items.append({
            "label": "Doublet removal policy",
            "value": "auto-applied (ATAC detector only; raw calls preserved)",
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
    items.append({
        "label": "Output location",
        "value": outputs,
        "reason": f"workflow_branch={branch}; run_manifest.json written alongside as handoff artifact.",
        "certainty": "certain",
    })

    # 8. Missing / uncertain
    missing = [k for k, v in ctx.get("fields", {}).items() if isinstance(v, dict) and v.get("status") == "missing"]
    low_conf = [k for k, v in ctx.get("fields", {}).items() if isinstance(v, dict) and v.get("confidence") == "low"]
    conflicts = ctx.get("conflicts", [])
    parts = []
    if missing:
        parts.append(f"missing: {', '.join(missing)}")
    if low_conf:
        parts.append(f"low-confidence: {', '.join(low_conf)}")
    if conflicts:
        parts.append(f"conflicts: {len(conflicts)}")
    items.append({
        "label": "Missing / uncertain info",
        "value": "; ".join(parts) if parts else "none",
        "reason": "Derived from P1 context field statuses and conflict list.",
        "certainty": "certain" if not parts else "needs confirmation",
    })

    return items


def render_summary_text(items: list[dict[str, Any]]) -> str:
    lines: list[str] = ["# Preprocessing plan review (concise)", ""]
    for it in items:
        tag = "✓" if it["certainty"] == "certain" else "?"
        lines.append(f"- **{it['label']}** [{tag} {it['certainty']}]")
        lines.append(f"  - value: `{it['value']}`")
        if it["reason"]:
            lines.append(f"  - reason: {it['reason']}")
    lines.append("")
    lines.append("_Full parameter details: `parameters.yaml` + `artifacts/p2_plan/preprocessing_plan.json`_")
    return "\n".join(lines)


def write_summary(run_dir: Path | str) -> Path:
    """Write the plan-review markdown directly to its canonical deliverable path."""
    from .run_paths import RunPaths
    items = build_summary(run_dir)
    out = RunPaths(Path(run_dir)).plan_review_md
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render_summary_text(items))
    return out
