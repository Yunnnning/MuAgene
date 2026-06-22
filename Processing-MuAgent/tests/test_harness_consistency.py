"""Harness consistency tripwires.

These guard the single-source-of-truth invariants introduced by the agent-harness
refactor. They are pure-Python and fast (no pipeline run). If one fails, a value
or contract has drifted from its canonical home — fix the source, not the test.

Stage 1 (this file's initial scope): every QC default lives once in
``executor.defaults.QC_DEFAULTS``. The plan assembler (which writes
``preprocessing_plan.json``) and ``executor.figures``' ``DEFAULT_*`` reference
constants (used by the pre-plan ``qc_explore`` preview) must both read from it, so
the plan, the stages, and the preview can never silently disagree.
"""
from __future__ import annotations

from executor import defaults
from executor import plan_assembler as pa


def test_plan_assembler_values_match_qc_defaults(tmp_path):
    """assemble_plan must emit exactly the centralised QC_DEFAULTS values+types.

    `paired` includes every QC-bearing stage. A literal sneaking back into
    plan_assembler (instead of reading QC_DEFAULTS) breaks this.
    """
    plan = pa.assemble_plan(tmp_path, workflow_branch="paired")
    stages = plan["stages"]
    for stage, params in defaults.QC_DEFAULTS.items():
        assert stage in stages, f"{stage} missing from assembled plan"
        for name, expected in params.items():
            got = stages[stage]["parameters"][name]["value"]
            assert got == expected, f"{stage}.{name}: plan={got!r} != defaults={expected!r}"
            # type matters: preprocessing_plan.json serialises int floors as `500`,
            # not `500.0` — a type drift would change the artifact byte-for-byte.
            assert type(got) is type(expected), (
                f"{stage}.{name}: type drift plan={type(got).__name__} "
                f"defaults={type(expected).__name__}")


def test_figures_default_constants_match_qc_defaults():
    """figures.DEFAULT_* (re-exported from defaults, consumed by qc_explore) must
    equal QC_DEFAULTS. Floors are exposed as float for marker geometry."""
    from executor import figures as F

    d = defaults.QC_DEFAULTS
    rna, atac = d["s1_rna_qc"], d["s2_atac_qc"]
    assert F.DEFAULT_TOTAL_COUNTS_K_MAD == rna["total_counts_k_mad"]
    assert F.DEFAULT_N_GENES_K_MAD == rna["n_genes_k_mad"]
    assert F.DEFAULT_PCT_MT_K == rna["pct_mt_k"]
    assert F.DEFAULT_PCT_MT_CEILING == rna["pct_mt_ceiling"]
    assert F.DEFAULT_PCT_MT_FLOOR == rna["pct_mt_floor"]
    assert F.DEFAULT_PCT_RIBO_MAX == rna["pct_ribo_max"]
    assert F.DEFAULT_MIN_COUNTS_FLOOR == float(rna["min_counts_floor"])
    assert F.DEFAULT_MIN_GENES_FLOOR == float(rna["min_genes_floor"])
    assert F.DEFAULT_N_FRAG_K_MAD == atac["n_fragments_k_mad"]
    assert F.DEFAULT_N_FRAG_FLOOR == float(atac["n_fragments_floor"])
    assert F.DEFAULT_TSS_MIN == atac["tss_enrichment_min"]
    assert F.DEFAULT_TSS_MAX == atac["tss_enrichment_max"]
    assert F.DEFAULT_NUC_MAX == atac["nucleosome_signal_max"]
