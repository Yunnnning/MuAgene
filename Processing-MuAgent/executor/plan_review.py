"""User-facing plan review: concise summary built from P1 context + S0 ingest + P2 plan."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import provenance


def _load_json(p: Path) -> dict:
    return json.loads(p.read_text()) if p.exists() else {}


# --- Marker-gene-check decision gate -------------------------------------
# Whenever ambient RNA correction is planned and no marker genes are set, the
# user must make an explicit choice before plan_review can be approved: provide
# genes, defer the check to QC review, or decline it. This turns the rendered
# "needs confirmation" item into an enforced gate so the check can never be
# silently skipped (e.g. by `submit --auto-approve`).

MARKER_GENE_DECISIONS = {"provided", "deferred_to_qc", "declined"}


def _planned_ambient_method(params_path: Path, plan: dict) -> str:
    """Effective S1a method: parameters.yaml (revise) wins over the frozen plan."""
    v = provenance.get_value(params_path, "s1a_ambient.method", None)
    if v:
        return str(v)
    p = plan.get("stages", {}).get("s1a_ambient", {}).get("parameters", {})
    m = p.get("method", {})
    return str(m.get("value", "auto")) if isinstance(m, dict) else "auto"


def marker_gene_decision_pending(run_dir: Path | str) -> bool:
    """True iff ambient correction is planned, no marker genes are set, AND no
    explicit decision (provided / deferred_to_qc / declined) has been recorded."""
    from .run_paths import RunPaths
    from .stages.s1a_ambient import resolve_marker_genes

    paths = RunPaths(Path(run_dir))
    plan = _load_json(paths.preprocessing_plan)
    if _planned_ambient_method(paths.parameters_yaml, plan).lower() in ("none", "skipped_empty"):
        return False
    plan_s1a = plan.get("stages", {}).get("s1a_ambient", {}).get("parameters", {})
    if resolve_marker_genes(paths.parameters_yaml, plan_s1a):
        return False  # genes provided — nothing pending
    decision = provenance.get_value(
        paths.parameters_yaml, "s1a_ambient.marker_genes_decision", None)
    return decision not in MARKER_GENE_DECISIONS


def record_marker_gene_decision(run_dir: Path | str, decision: str) -> None:
    """Persist an explicit user decision to defer or decline the marker check."""
    if decision not in ("deferred_to_qc", "declined"):
        raise ValueError(
            f"marker gene decision must be 'deferred_to_qc' or 'declined', got {decision!r}")
    from .run_paths import RunPaths

    paths = RunPaths(Path(run_dir))
    rationale = (
        "User chose to check marker genes at QC review instead of plan review."
        if decision == "deferred_to_qc"
        else "User declined the marker gene expression check for ambient correction."
    )
    provenance.set_param(
        paths.parameters_yaml, "s1a_ambient.marker_genes_decision", decision,
        source="user", confidence="high", rationale=rationale,
    )


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
    ingest = _load_json(paths.validation_report)
    plan = _load_json(paths.preprocessing_plan)

    def field(name: str) -> dict[str, Any]:
        return ctx.get("fields", {}).get(name, {})

    pairing = ingest.get("pairing", {})
    overlap_raw = pairing.get("overlap")
    return {
        "organism":             field("organism").get("value", ""),
        "tissue":               field("tissue").get("value", ""),
        "assay_type":           field("assay_type").get("value", ""),
        "sample_type":          field("sample_type").get("value", ""),
        "genome_build":         field("genome_build").get("value", ""),
        "workflow_branch":      plan.get("workflow_branch", ingest.get("workflow_branch", "")),
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
    from .plan_assembler import overlay_plan
    ctx = _load_json(paths.artifact("p1_context", "context_extraction.json"))
    ingest = _load_json(paths.validation_report)
    # Effective plan: frozen plan overlaid with any parameters.yaml `revise`, so
    # the review shows what stages will actually apply.
    plan = overlay_plan(_load_json(paths.preprocessing_plan), paths.parameters_yaml)
    pairing = ingest.get("pairing", {})

    def field(name: str) -> dict[str, Any]:
        return ctx.get("fields", {}).get(name, {})

    def param(stage: str, key: str) -> dict[str, Any]:
        return plan.get("stages", {}).get(stage, {}).get("parameters", {}).get(key, {})

    items: list[dict[str, Any]] = []

    # 1. Execution mode + HPC settings (from configure-execution / hpc.env)
    items.extend(_execution_review_items(run_dir, plan))

    # 2. Dataset type / modality
    modality = field("modality_type")
    _pdeclared = pairing.get("declared", "")
    _pcommitted = pairing.get("committed", "")
    _pconfidence = pairing.get("confidence", "")
    _branch_conflict = bool(_pdeclared and _pcommitted and _pdeclared != _pcommitted)
    if modality.get("status") == "explicit":
        _dt_certainty = "certain"
        _dt_value = modality.get("value", "unknown")
        _dt_reason = modality.get("rationale", "")
    elif _branch_conflict:
        _dt_certainty = "conflict"
        _dt_value = modality.get("value", "unknown")
        _dt_reason = (
            f"Declared as '{_pdeclared}' but S0 investigation committed '{_pcommitted}'. "
            + modality.get("rationale", "")
        ).strip()
    elif _pconfidence == "high":
        _dt_certainty = "certain"
        _dt_value = plan.get("workflow_branch", modality.get("value", "unknown"))
        _dt_reason = _pairing_review_reason(pairing)
    else:
        _dt_certainty = "needs confirmation"
        _dt_value = modality.get("value", "unknown")
        _dt_reason = modality.get("rationale", "")
    items.append({
        "label": "Detected dataset type",
        "value": _dt_value,
        "reason": _dt_reason,
        "certainty": _dt_certainty,
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
    items.append({
        "label": "Pairing (RNA ↔ ATAC)",
        "value": f"{pairing.get('status', 'unknown')} (overlap={pairing.get('overlap', 0):.2%})",
        "reason": _pairing_review_reason(pairing),
        "certainty": (
            "conflict" if _branch_conflict
            else "certain" if pairing.get("confidence") == "high"
            else "needs confirmation"
        ),
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

    # 5. Marker gene expression check — explicit user confirmation required
    from .stages.s1a_ambient import resolve_marker_genes

    plan_s1a_params = plan.get("stages", {}).get("s1a_ambient", {}).get("parameters", {})
    mg_list = resolve_marker_genes(paths.parameters_yaml, plan_s1a_params)
    mg_display = ", ".join(mg_list) if mg_list else "not set"
    items.append({
        "label": "Marker gene expression check",
        "value": mg_display,
        "reason": (
            "Would you like to visualise how 5–10 marker genes distribute across "
            "cell clusters before and after Ambient RNA Correction? "
            "If a marker gene appears at low levels ubiquitously across cells that "
            "shouldn't express it, this is a sign of ambient RNA contamination. "
            "After correction, expression should be clearer and more restricted to "
            "the expected populations. "
            "If yes, provide gene symbols (e.g. CD3E, CD20, EPCAM). "
            "If no, leave as not set — the check is skipped."
        ),
        "certainty": "needs confirmation",
    })

    # 6. Key QC strategy
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
        "reason": (
            "Review the QC threshold histograms in the appendix. Default MAD thresholds "
            "are shown — any RNA or ATAC metric can be adjusted or skipped entirely with "
            "`revise` before approving. Confirm defaults are acceptable, or tell the agent "
            "which thresholds to change."
        ),
        "certainty": "needs confirmation",
    })

    # 5. Doublet policy — reconciliation is only meaningful when both detectors run.
    # Single-modality branches have one detector; the S3 item is informational, not
    # a checkpoint gate.
    branch = plan.get("workflow_branch", "paired")
    policy = param("s3_doublets", "removal_policy_recommendation")

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
            "certainty": "needs confirmation",
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
            "certainty": "needs confirmation",
        })
    elif branch == "rna_only":
        rna_thr_note = ""
        if rna_doub_thr.get("value") is not None:
            rna_thr_note = f"; score threshold={rna_doub_thr.get('value')}"
        items.append({
            "label": "Doublet removal policy",
            "value": f"Scrublet only (fixed score threshold{rna_thr_note})",
            "reason": "single-detector branch — no reconciliation to confirm.",
            "certainty": "needs confirmation",
        })
    elif branch == "atac_only":
        atac_thr_note = ""
        if atac_doub_thr.get("value") is not None:
            atac_thr_note = f"; score threshold={atac_doub_thr.get('value')}"
        items.append({
            "label": "Doublet removal policy",
            "value": f"SnapATAC2 scrublet only (fixed score threshold{atac_thr_note})",
            "reason": "single-detector branch — no reconciliation to confirm.",
            "certainty": "needs confirmation",
        })

    # 6. Clustering
    rna_res = param("s7_clustering", "rna_resolution")
    atac_res = param("s7_clustering", "atac_resolution")
    items.append({
        "label": "Clustering strategy",
        "value": f"Leiden at fixed resolutions (RNA={rna_res.get('value', '?')}, ATAC={atac_res.get('value', '?')})",
        "reason": "Fixed per-modality defaults; clustering runs automatically with no resolution checkpoint.",
        "certainty": "needs confirmation",
    })

    # 7. Output location
    branch = plan.get("workflow_branch", "unknown")
    outputs = {
        "paired": str(paths.processed_h5mu),
        "separate": f"{paths.rna_processed_h5ad} + {paths.atac_processed_h5ad}",
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
        "certainty": "conflict" if conflicts else ("certain" if not parts else "needs confirmation"),
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
            "How heavy stages (S0 onward) are dispatched: local machine vs SLURM.",
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

    # SLURM — surface scheduler settings for review
    parts: list[str] = []
    missing: list[str] = []
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
        tag = {"certain": "✓", "conflict": "⚠"}.get(it["certainty"], "?")
        lines.append(f"- **{it['label']}** [{tag} {it['certainty']}]")
        lines.append(f"  - value: `{it['value']}`")
        for detail_line in it.get("detail", []):
            lines.append(f"  {detail_line}")
        if it.get("reason"):
            lines.append(f"  - reason: {it['reason']}")
    return "\n".join(lines)


def plan_review_heading(run_dir: Path | str) -> str:
    """User-facing plan-review title including the run name."""
    return f"Preprocessing plan review — {Path(run_dir).name}"


def _persist_intro(paths, intro: str) -> None:
    """Persist the agent-authored intro paragraph next to the plan.

    Lets every later render (the propose rule, HTML regeneration, resume) reuse
    the same intro instead of dropping it — it is otherwise only a transient
    ``--intro`` CLI argument.
    """
    text = (intro or "").strip()
    if not text:
        return
    paths.plan_intro.parent.mkdir(parents=True, exist_ok=True)
    paths.plan_intro.write_text(text + "\n")


def _load_persisted_intro(paths) -> str | None:
    """Return the persisted intro paragraph, or None if none was stored."""
    p = paths.plan_intro
    if p.exists():
        text = p.read_text().strip()
        return text or None
    return None


def render_merged_markdown(run_dir: Path | str, intro: str | None = None) -> str:
    """Full plan-review deliverable: concise summary + parameter appendix.

    When ``intro`` is None, falls back to the persisted intro paragraph (written
    by ``executor plan-review --intro``) so propose re-renders keep it.
    """
    from .plan_assembler import overlay_plan, render_plan_appendix
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    paths = RunPaths(run_dir)
    if intro is None:
        intro = _load_persisted_intro(paths)
    items = build_summary(run_dir)
    # Effective plan: frozen plan overlaid with any parameters.yaml `revise`, so
    # the appendix matches what stages will apply (the QC blocks reflect it too).
    plan = overlay_plan(_load_json(paths.preprocessing_plan), paths.parameters_yaml)
    parts = [f"# {plan_review_heading(run_dir)}", ""]
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
        from . import qc_explore
        # The Snakemake plan_review rule builds qc_explore.json first; when this is
        # called directly (`executor plan-review`) compute it once, best-effort.
        if not paths.artifact("qc_explore", "qc_explore.json").exists():
            try:
                qc_explore.run(run_dir)
            except Exception:
                pass
        qc_blocks = qc_explore.render_appendix_blocks(run_dir)
        parts.append(render_plan_appendix(plan, qc_blocks).rstrip())
        parts.append("")
    else:
        parts.append("_Appendix unavailable — `preprocessing_plan.json` not found._")
        parts.append("")
    return "\n".join(parts)


def _remove_legacy_plan_review_paths(paths) -> None:
    """Drop pre-run-scoped filenames (plan_review.md / plan_summary.html)."""
    for legacy in (
        paths.deliv_plan / "plan_review.md",
        paths.deliv_plan / "plan_summary.html",
    ):
        if legacy.exists():
            legacy.unlink()


def write_summary(run_dir: Path | str, intro: str | None = None) -> Path:
    """Write merged plan-review markdown (summary + parameter appendix)."""
    from .run_paths import RunPaths
    run_dir = Path(run_dir)
    paths = RunPaths(run_dir)
    if intro is not None:
        _persist_intro(paths, intro)
    merged = render_merged_markdown(run_dir, intro=intro)
    paths.plan_review_md.parent.mkdir(parents=True, exist_ok=True)
    _remove_legacy_plan_review_paths(paths)
    paths.plan_review_md.write_text(merged)
    return paths.plan_review_md


def render_plan_summary_html(run_dir: Path | str, intro: str | None = None) -> str:
    """Self-contained HTML rendering of the plan review.

    Same content as ``plan_review.md`` but figures are embedded as base64 data
    URIs so the single file is viewable on HPC nodes without figure-preview
    support. Reuses the markdown→HTML helpers from ``qc_summary``.
    """
    from .qc_summary import _embed_html_images, _markdown_to_html, _qc_html_styles
    from .run_paths import RunPaths
    paths = RunPaths(Path(run_dir))
    heading = plan_review_heading(run_dir)
    md = render_merged_markdown(run_dir, intro=intro)
    body = _markdown_to_html(md)
    # Image src in the markdown is relative to plan_review.md's directory.
    body = _embed_html_images(body, paths.plan_review_md.parent)
    return (
        "<!DOCTYPE html>\n<html lang=\"en\"><head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        f"<title>{heading}</title>\n"
        f"<style>\n{_qc_html_styles()}</style>\n</head>\n<body>\n"
        f"{body}\n</body></html>\n"
    )


def write_plan_summary_html(run_dir: Path | str, intro: str | None = None) -> Path:
    """Write the self-contained HTML plan review next to plan_review.md."""
    from .run_paths import RunPaths
    paths = RunPaths(Path(run_dir))
    if intro is not None:
        _persist_intro(paths, intro)
    html = render_plan_summary_html(run_dir, intro=intro)
    paths.plan_summary_html.parent.mkdir(parents=True, exist_ok=True)
    _remove_legacy_plan_review_paths(paths)
    paths.plan_summary_html.write_text(html)
    return paths.plan_summary_html
