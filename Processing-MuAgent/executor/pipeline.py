"""Pipeline topology — the single source of truth for stage identity and order.

Every other module derives its view of "which stages exist, in what order, and which
run on a given workflow branch" from the constants here. Before this module existed the
same lists lived (and drifted) in three places: cli.STAGES, plan_assembler._STAGES_BY_BRANCH,
and specs._BRANCH_STAGES — a drift that left s6_neighbors/s7_clustering omitted from the
atac_only specs even though the Snakemake DAG runs them on every branch. Keep this module
authoritative: it must agree with the actual Snakemake DAG (workflow/rules/*.smk), and a
harness-consistency test asserts the consumers stay in sync.

Branch membership reflects what the DAG executes:
  - paired / unpaired : the full modality pipeline (RNA + ATAC).
  - rna_only          : drops the ATAC-only stages (s2_atac_qc, s5_atac_spectral).
  - atac_only         : drops the RNA-only stages (s1a_ambient, s1_rna_qc, s4_rna_norm).
s3_doublets, s6_neighbors, s7_clustering, s8_umap run on ALL branches — s6/s7/s8 have no
branch guard in the DAG (s6 keys off the S5 spectral marker on atac_only), so they cluster
whichever modality is present.
"""
from __future__ import annotations

from .provenance import branch_has_atac  # re-export: branch → has ATAC modality?

__all__ = [
    "STAGES",
    "PLANNING_STAGES",
    "PIPELINE_STAGE_ORDER",
    "HUMAN_CHECKPOINTS",
    "QC_HANDOFF",
    "AUTOMATED_STAGES",
    "STAGES_BY_BRANCH",
    "STAGE_ALIASES",
    "STAGE_DISPLAY",
    "branch_has_atac",
    "stages_for_branch",
    "canonical_stage",
    "display_stage",
]

# Planning compute that precedes the modality pipeline. s0_ingest is the merged
# planning stage (load + validate + assemble plan + QC exploration); the former
# standalone p2_plan stage no longer exists.
PLANNING_STAGES: tuple[str, ...] = ("p1_context", "s0_ingest", "plan_review")

# The modality pipeline in execution order (S1a before S1, etc.). This is the set of
# stages that appear under plan["stages"] and is the basis for per-branch membership.
PIPELINE_STAGE_ORDER: tuple[str, ...] = (
    "s1a_ambient",
    "s1_rna_qc",
    "s2_atac_qc",
    "s3_doublets",
    "s4_rna_norm",
    "s5_atac_spectral",
    "s6_neighbors",
    "s7_clustering",
    "s8_umap",
)

# The two human gates (Snakemake sentinels).
HUMAN_CHECKPOINTS: tuple[str, ...] = ("plan_review", "post_qc_review")

# The post-QC Integration bundle — a cluster job that runs on every branch but is not
# part of plan["stages"] (it has no per-branch membership entry of its own).
QC_HANDOFF: str = "qc_handoff"

# Every stage id the CLI accepts, in workflow order.
STAGES: tuple[str, ...] = (
    *PLANNING_STAGES,
    "s1a_ambient", "s1_rna_qc", "s2_atac_qc", "s3_doublets",
    "post_qc_review", QC_HANDOFF,
    "s4_rna_norm", "s5_atac_spectral", "s6_neighbors", "s7_clustering", "s8_umap",
)

# Stages that advance without a human gate.
AUTOMATED_STAGES: tuple[str, ...] = tuple(s for s in STAGES if s not in HUMAN_CHECKPOINTS)

# Per-branch membership for the modality pipeline (the stages that appear under
# plan["stages"] and get per-stage specs). MUST match the Snakemake DAG.
STAGES_BY_BRANCH: dict[str, set[str]] = {
    "paired":    set(PIPELINE_STAGE_ORDER),
    "unpaired":  set(PIPELINE_STAGE_ORDER),
    "rna_only":  {"s1a_ambient", "s1_rna_qc", "s3_doublets", "s4_rna_norm",
                  "s6_neighbors", "s7_clustering", "s8_umap"},
    "atac_only": {"s2_atac_qc", "s3_doublets", "s5_atac_spectral",
                  "s6_neighbors", "s7_clustering", "s8_umap"},
}

# Stage aliases: the user-/agent-facing name `qc_review` maps to the internal stage
# `post_qc_review`; display maps back for human-readable echoes.
STAGE_ALIASES: dict[str, str] = {
    "qc_review": "post_qc_review",
}
STAGE_DISPLAY: dict[str, str] = {
    "post_qc_review": "qc_review",
}


def stages_for_branch(branch: str) -> set[str]:
    """Return the modality stage ids that run for `branch` (KeyError → ValueError)."""
    try:
        return STAGES_BY_BRANCH[branch]
    except KeyError as exc:
        raise ValueError(
            f"Unknown workflow_branch={branch!r}; expected one of {sorted(STAGES_BY_BRANCH)}."
        ) from exc


def canonical_stage(stage: str) -> str:
    """Map a user-facing stage alias to its internal stage id (passthrough if none)."""
    return STAGE_ALIASES.get(stage, stage)


def display_stage(stage: str) -> str:
    """Map an internal stage id to its user-facing display name (passthrough if none)."""
    return STAGE_DISPLAY.get(stage, stage)
