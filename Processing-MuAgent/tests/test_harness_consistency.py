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

import json
import pathlib

import pytest

from executor import defaults
from executor import plan_assembler as pa


def _contracts_dir() -> pathlib.Path:
    # Processing-MuAgent/tests/<this> -> parents[2] == MuAgene repo root.
    return pathlib.Path(__file__).resolve().parents[2] / "contracts"


def _repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[2]


def test_root_agent_is_the_two_agent_entry_point():
    root_agent = (_repo_root() / "AGENT.md").read_text()
    assert "Processing-MuAgent/AGENT.md" in root_agent
    assert "Execution-MuAgent/AGENT.md" in root_agent
    assert "muagene.agents.yaml" in root_agent
    assert "## Global safeguards" not in root_agent
    import yaml
    registry = yaml.safe_load((_repo_root() / "muagene.agents.yaml").read_text())
    assert registry["root_manifest"] == "AGENT.md"


def test_component_instruction_hierarchies_start_at_root_agent():
    for component in ("Processing-MuAgent", "Execution-MuAgent"):
        manifest = (_repo_root() / component / "AGENT.md").read_text()
        router = (_repo_root() / component / "agent" / "skills" / "index.md").read_text()
        loading_order = router.split("## Loading order", 1)[1].split("\n## ", 1)[0]
        assert "Policy + entry point" not in manifest
        assert loading_order.index("../../../AGENT.md") < loading_order.index("../../AGENT.md")


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


# --- Stage 2: contracts/post_qc_manifest.schema.json ---

def _representative_manifest() -> dict:
    """A manifest with exactly the keys/types qc_handoff.run() emits."""
    from executor import HANDOFF_CONTRACT_VERSION
    return {
        "schema": "muagene.post_qc_handoff/1",
        "handoff_contract_version": HANDOFF_CONTRACT_VERSION,
        "sample_run_dir": "/runs/sampleA",
        "modality_branch": "paired",
        "genome_assembly": "GRCh38",
        "post_qc_h5mu": "deliverables/qc/post_qc_sampleA.h5mu",
        "atac": {
            "peaks_bed": "deliverables/qc/peaks_sampleA.bed",
            "peaks_source": "macs3",
            "fragments_prepared": "internal/artifacts/s2_atac_qc/atac_fragments_cbf.tsv.gz",
            "add_chr_prefix": True,
            "frag_chrom_convention": "ucsc",
        },
        "n_cells": {"rna": 100, "atac": 90, "joint": 95},
        "parameters_ref": "internal/parameters.yaml",
        "tool_versions": {"scanpy": "1.10.0"},
    }


def test_post_qc_manifest_schema_is_wellformed():
    schema = json.loads((_contracts_dir() / "post_qc_manifest.schema.json").read_text())
    assert schema["$id"] == "muagene.post_qc_handoff/1"
    assert "schema" in schema["required"]
    assert schema["properties"]["modality_branch"]["enum"] == [
        "paired", "unpaired", "rna_only", "atac_only"]
    # Every schema-required top-level key is one the emitter actually writes.
    assert set(schema["required"]) <= set(_representative_manifest())


def test_post_qc_manifest_representative_validates_against_schema():
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads((_contracts_dir() / "post_qc_manifest.schema.json").read_text())
    jsonschema.validate(_representative_manifest(), schema)  # valid -> no raise
    bad = _representative_manifest()
    bad["modality_branch"] = "bogus"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(bad, schema)


def test_run_manifest_contract_version_matches_emitter():
    import yaml
    from executor import HANDOFF_CONTRACT_VERSION

    registry = yaml.safe_load(
        (pathlib.Path(__file__).resolve().parents[2] / "muagene.agents.yaml").read_text()
    )
    assert registry["handoff_contracts"]["run_manifest"] == HANDOFF_CONTRACT_VERSION


# --- Stage 5: every CLI command has a tool contract ---

def test_every_executor_command_is_documented():
    import re
    from executor.cli import main
    tools = (pathlib.Path(__file__).resolve().parents[1] / "agent" / "tools.md").read_text()
    documented = set(re.findall(r"^### executor ([a-z0-9-]+)$", tools, flags=re.MULTILINE))
    live = set(main.commands)
    assert documented == live, (
        f"agent/tools.md mismatch: missing={sorted(live - documented)}, "
        f"stale={sorted(documented - live)}"
    )


# --- Stage 6: revise --dry-run previews without mutating ---

def test_qc_downstream_targets_is_nonmutating(tmp_path):
    """The preview helper behind `revise --dry-run` lists the would-delete artifacts
    but deletes nothing — the safeguard against a destructive revise."""
    from executor import cli
    from executor.run_paths import RunPaths
    rp = RunPaths(tmp_path)
    f = rp.artifact("s3_doublets", "calls.parquet")
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text("x")
    targets = cli._qc_downstream_targets(tmp_path, "s1_rna_qc")
    assert f in targets       # an s1 revise would invalidate the downstream s3 artifact
    assert f.exists()         # but computing the preview deletes nothing


def test_revise_has_dry_run_flag():
    from executor.cli import main
    params = {p.name for p in main.commands["revise"].params}
    assert "dry_run" in params, f"revise is missing the --dry-run flag; has {sorted(params)}"


# --- Stage 8: stage-based skill router + per-skill frontmatter contracts ---

def _skills_dir() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1] / "agent" / "skills"


# The contract keys every skill's YAML frontmatter must carry (index.md is the router,
# not a skill, and is exempt). Mirrors the schema documented in agent/skills/index.md.
_REQUIRED_FRONTMATTER = {
    "name", "domain", "purpose", "activation", "inputs", "outputs",
    "calls_tools", "reads_contracts", "writes_state", "handoff",
}

_ORDERED_SKILLS = {
    "00_entry_declare.md": "entry_declare",
    "10_inputs_intake.md": "inputs_intake",
    "20_plan_confirm.md": "plan_confirm",
    "30_run_execution.md": "run_execution",
    "40_qc_review_and_revise.md": "qc_review_and_revise",
    "50_downstream_dimred_clustering.md": "downstream_dimred_clustering",
    "60_completion_handoff.md": "completion_handoff",
    "80_hpc_monitoring.md": "hpc_monitoring",
    "90_troubleshooting.md": "troubleshooting",
}


def test_processing_skill_filenames_encode_router_order():
    actual = sorted(md.name for md in _skills_dir().glob("*.md") if md.name != "index.md")
    assert actual == list(_ORDERED_SKILLS)


def test_processing_skill_frontmatter_preserves_semantic_ids():
    import yaml
    for filename, semantic_id in _ORDERED_SKILLS.items():
        text = (_skills_dir() / filename).read_text()
        frontmatter = yaml.safe_load(text.split("---\n", 2)[1])
        assert isinstance(frontmatter, dict), filename
        assert _REQUIRED_FRONTMATTER <= set(frontmatter), filename
        assert frontmatter["name"] == semantic_id, filename


def test_qc_review_skill_is_compact_and_preserves_safety_actions():
    text = (_skills_dir() / "40_qc_review_and_revise.md").read_text()
    required = [
        "--dry-run",
        "Confirm with the user",
        "do not hand-delete",
        "Do not use `--auto-approve`",
        "relay it **verbatim**",
        "never pick genes",
        "executor marker-gene-check --config $CFG",
        "executor approve post_qc_review --config $CFG",
        "--target s3_doublets_execute",
        "executor propose post_qc_review --config $CFG",
        "--target qc_handoff",
        "post_qc_<run>.h5mu",
        "post_qc_manifest.json",
        "explicit user approval",
        "total_counts_min_override",
        "pct_mt_k=999",
    ]
    missing = [phrase for phrase in required if phrase not in text]
    assert not missing, f"QC skill lost safety actions: {missing}"
    assert len(text.splitlines()) <= 220
    assert "executor approve s3_doublets" not in text
    assert "executor approve <stage>" not in text


def test_no_markdown_references_unprefixed_processing_skill_filenames():
    old_names = [filename.split("_", 1)[1] for filename in _ORDERED_SKILLS]
    for markdown in _repo_root().rglob("*.md"):
        text = markdown.read_text()
        stale = [name for name in old_names if f"({name})" in text or f"skills/{name}" in text]
        assert not stale, f"{markdown.relative_to(_repo_root())}: stale skill links {stale}"


def test_every_skill_has_required_frontmatter():
    """Each stage skill opens with a frontmatter contract carrying every required key.
    A new skill that forgets purpose/activation/handoff (the routing+contract fields)
    fails here."""
    for md in sorted(_skills_dir().glob("*.md")):
        if md.name == "index.md":
            continue
        text = md.read_text()
        assert text.startswith("---\n"), f"{md.name}: missing YAML frontmatter block"
        front = text.split("---\n", 2)[1]
        top_keys = {
            line.split(":", 1)[0].strip()
            for line in front.splitlines()
            if ":" in line and not line.startswith((" ", "\t"))
        }
        missing = _REQUIRED_FRONTMATTER - top_keys
        assert not missing, f"{md.name}: frontmatter missing {sorted(missing)}"


def test_skill_cross_links_resolve():
    """Every relative .md link inside a skill points to a file that exists — guards
    against dangling pointers after the workflow.md -> stage-skills split."""
    import re
    link_re = re.compile(r"\]\(([^)]+\.md)\)")
    for md in sorted(_skills_dir().glob("*.md")):
        for target in link_re.findall(md.read_text()):
            if target.startswith("http"):
                continue
            resolved = (md.parent / target).resolve()
            assert resolved.exists(), f"{md.name}: dangling link -> {target}"


def test_router_lists_every_skill():
    """index.md (the router) names every sibling skill; the router and the skill set
    cannot drift apart."""
    index = (_skills_dir() / "index.md").read_text()
    for md in sorted(_skills_dir().glob("*.md")):
        if md.name == "index.md":
            continue
        assert md.name in index, f"router index.md does not list {md.name}"


def test_no_skill_references_deleted_workflow_md():
    """workflow.md was dissolved into stage skills; nothing should reference it again."""
    for md in sorted(_skills_dir().glob("*.md")):
        assert "workflow.md" not in md.read_text(), f"{md.name} still references workflow.md"


# --- Stage 9: per-stage specs mirror the pipeline SSOT and the DAG's durable markers ---

def test_specs_membership_matches_pipeline_ssot():
    """The stages specs writes for a branch (minus the always-present s0_ingest planning
    compute + qc_handoff bundle) must equal pipeline.stages_for_branch — so specs can never
    re-introduce the atac_only s6/s7 drift this refactor fixed."""
    from executor import pipeline, specs
    for branch in pipeline.STAGES_BY_BRANCH:
        spec_modality = set(specs._spec_stages(branch)) - {"s0_ingest", "qc_handoff"}
        assert spec_modality == pipeline.stages_for_branch(branch), branch


def test_atac_only_specs_include_s6_s7_and_exclude_rna_stages():
    """Ground truth: the Snakemake DAG runs s6_neighbors/s7_clustering on every branch."""
    from executor import specs
    stages = set(specs._spec_stages("atac_only"))
    assert {"s6_neighbors", "s7_clustering"} <= stages
    assert not ({"s1a_ambient", "s1_rna_qc", "s4_rna_norm"} & stages)


def test_stage_spec_outputs_are_durable_markers():
    """Each per-stage spec advertises its durable stage-done marker as the declared
    output — what the monitor can verify and what survives cleanup — not a deletable
    working h5ad. The marker set is stage_progress.EXECUTE_MARKERS (the DAG's edge keys)."""
    from executor import specs
    from executor.stage_progress import EXECUTE_MARKERS
    for stage, marker in EXECUTE_MARKERS.items():
        _inputs, outputs = specs._stage_io_for_branch(stage, "paired", "/runs/x")
        vals = list(outputs.values())
        assert any(v.endswith(f"{stage}/{marker}") for v in vals), (stage, marker, vals)


def test_s6_inputs_are_branch_aware():
    """s6_neighbors keys off the S4 RNA marker on RNA branches and the S5 spectral marker
    on ATAC branches — exactly mirroring s6_neighbors.smk's _s6_inputs."""
    from executor import specs
    rna_in, _ = specs._stage_io_for_branch("s6_neighbors", "rna_only", "/runs/x")
    atac_in, _ = specs._stage_io_for_branch("s6_neighbors", "atac_only", "/runs/x")
    assert any("s4_rna_norm/norm_summary.json" in v for v in rna_in.values())
    assert not any("s5_atac_spectral" in v for v in rna_in.values())
    assert any("s5_atac_spectral/spectral_summary.json" in v for v in atac_in.values())
    assert not any("s4_rna_norm" in v for v in atac_in.values())
