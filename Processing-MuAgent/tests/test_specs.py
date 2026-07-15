"""Unit tests for executor.specs — per-stage job specs handed to Execution-MuAgent.

Focus: the spec set + I/O follow the pipeline SSOT and the DAG's durable markers,
including the atac_only fix (s6/s7 now specced, keyed off the S5 spectral marker).
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import yaml

from executor import specs


def _written_stage_names(branch: str) -> set[str]:
    with tempfile.TemporaryDirectory() as tmp:
        written = specs.write_stage_specs(tmp, branch)
        return {p.stem for p in written}


def test_atac_only_writes_s6_s7_specs():
    names = _written_stage_names("atac_only")
    assert {"s6_neighbors", "s7_clustering"} <= names
    # RNA-only stages get no spec on atac_only.
    assert not ({"s1a_ambient", "s1_rna_qc", "s4_rna_norm"} & names)
    # The ATAC stages and the always-present brackets are present.
    assert {"s0_ingest", "s2_atac_qc", "s3_doublets", "s5_atac_spectral",
            "s8_umap", "qc_handoff"} <= names


def test_rna_only_writes_no_atac_stage_specs():
    names = _written_stage_names("rna_only")
    assert not ({"s2_atac_qc", "s5_atac_spectral"} & names)
    assert {"s6_neighbors", "s7_clustering"} <= names


def test_s6_spec_input_on_atac_only_is_spectral_marker():
    spec = _build("s6_neighbors", "atac_only")
    inputs = set(spec["inputs"].values())
    assert any(v.endswith("s5_atac_spectral/spectral_summary.json") for v in inputs)
    assert not any("s4_rna_norm" in v for v in inputs)
    # Output is the durable neighbors marker, not the deletable rna_neighbors.h5ad.
    assert any(v.endswith("s6_neighbors/neighbors_summary.json") for v in spec["outputs"].values())
    assert not any(v.endswith("rna_neighbors.h5ad") for v in spec["outputs"].values())


def test_s7_spec_io_uses_durable_markers():
    spec = _build("s7_clustering", "paired")
    assert any(v.endswith("s6_neighbors/neighbors_summary.json") for v in spec["inputs"].values())
    assert any(v.endswith("s7_clustering/clustering_summary.json") for v in spec["outputs"].values())


def test_qc_handoff_peaks_only_on_atac_branches():
    paired = _build("qc_handoff", "paired")
    rna = _build("qc_handoff", "rna_only")
    assert any(v.endswith(".bed") for v in paired["outputs"].values())
    assert not any(v.endswith(".bed") for v in rna["outputs"].values())


def _build(stage: str, branch: str) -> dict:
    """Write specs for a branch and return the parsed spec for one stage."""
    with tempfile.TemporaryDirectory() as tmp:
        specs.write_stage_specs(tmp, branch)
        path = Path(tmp) / "internal" / "stage_meta" / f"{stage}.yaml"
        return yaml.safe_load(path.read_text())
