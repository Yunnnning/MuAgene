"""Tests for the pre-plan QC exploration: non-exclusive counting, shared
threshold derivation, and plan-review table/figure embedding."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless figure rendering for the explore tests

import numpy as np
import pandas as pd

from executor import io as _io
from executor import qc_explore
from executor.methods.qc_filter_stats import marginal_removals
from executor.methods import qc_thresholds as qct
from executor.plan_assembler import render_plan_appendix
from executor.run_paths import RunPaths


class MarginalRemovalsTests(unittest.TestCase):
    def test_counts_are_non_exclusive(self):
        masks = {
            "a": np.array([True, False, False, True]),   # fails idx 1, 2
            "b": np.array([True, True, False, False]),   # fails idx 2, 3
        }
        rm = marginal_removals(masks)
        self.assertEqual(rm["a"], 2)            # every cell failing a
        self.assertEqual(rm["b"], 2)            # every cell failing b
        self.assertEqual(rm["multiple_metrics"], 1)  # idx 2 fails both
        self.assertEqual(rm["total_removed"], 3)     # union {1,2,3}

    def test_order_independent(self):
        a = np.array([True, False, False, True])
        b = np.array([True, True, False, False])
        self.assertEqual(
            marginal_removals({"a": a, "b": b}),
            marginal_removals({"b": b, "a": a}),
        )


class RnaThresholdTests(unittest.TestCase):
    def test_floors_clamp_lower_bounds(self):
        import pandas as pd
        rng = np.random.default_rng(0)
        n = 500
        obs = pd.DataFrame({
            "total_counts": rng.integers(50, 5000, n).astype(float),
            "n_genes_by_counts": rng.integers(20, 2000, n).astype(float),
            "pct_counts_mt": rng.uniform(0, 8, n),
            "pct_counts_ribo": rng.uniform(0, 30, n),
        })
        th = qct.rna_thresholds(
            obs, k_mad=5.0, pct_mt_k=3.0, pct_mt_ceiling=20.0, pct_mt_floor=5.0,
            min_counts_floor=500, min_genes_floor=200,
        )
        self.assertGreaterEqual(th["total_counts_min"], 500)
        self.assertGreaterEqual(th["n_genes_min"], 200)
        self.assertGreaterEqual(th["pct_counts_mt_max"], 5.0)
        self.assertIn("pct_counts_mt_mad_raw", th)
        self.assertLessEqual(th["pct_counts_mt_mad_raw"], th["pct_counts_mt_max"])

        # Pristine mito profile: MAD bound sits below the floor.
        pristine = pd.DataFrame({
            "total_counts": np.full(200, 2000.0),
            "n_genes_by_counts": np.full(200, 1500.0),
            "pct_counts_mt": rng.uniform(0, 3, 200),
        })
        th_p = qct.rna_thresholds(
            pristine, k_mad=5.0, pct_mt_k=3.0, pct_mt_ceiling=20.0, pct_mt_floor=5.0,
            min_counts_floor=500, min_genes_floor=200,
        )
        self.assertLess(th_p["pct_counts_mt_mad_raw"], th_p["pct_counts_mt_max"])
        masks = qct.rna_pass_masks(obs, th, pct_ribo_max=50.0)
        self.assertEqual(set(masks), {"total_counts", "n_genes",
                                       "pct_counts_mt", "pct_counts_ribo"})


class PctMtPanelRefsTests(unittest.TestCase):
    def test_shows_mad_when_floor_clamps_applied(self):
        th = {"pct_counts_mt_max": 5.0, "pct_counts_mt_mad_raw": 3.4}
        refs = qc_explore._pct_mt_panel_refs(th)
        labels = [r[1] for r in refs]
        self.assertIn("3.4% MAD", labels)
        self.assertIn("5%", labels)
        self.assertIn("10%", labels)

    def test_omits_mad_when_it_matches_applied(self):
        th = {"pct_counts_mt_max": 8.2, "pct_counts_mt_mad_raw": 8.2}
        refs = qc_explore._pct_mt_panel_refs(th)
        labels = [r[1] for r in refs]
        self.assertNotIn("8.2% MAD", labels)
        self.assertEqual(labels, ["5%", "10%"])


class AppendixBlockTests(unittest.TestCase):
    def _write_explore(self, run_dir: Path) -> RunPaths:
        paths = RunPaths(run_dir)
        art = paths.stage_dir("qc_explore")
        art.mkdir(parents=True, exist_ok=True)
        paths.deliv_figures.mkdir(parents=True, exist_ok=True)
        paths.plan_review_md.parent.mkdir(parents=True, exist_ok=True)
        (paths.deliv_figures / f"{qc_explore.S1_FIGURE_STEM}.png").write_bytes(b"x")
        data = {
            "s1_rna_qc": {
                "thresholds": {
                    "total_counts_min": 500.0, "total_counts_max": 40000.0,
                    "n_genes_min": 200.0, "n_genes_max": 8000.0,
                    "pct_counts_mt_max": 12.0, "pct_counts_ribo_max": 50.0,
                },
                "cells_removed": {
                    "total_counts": 100, "n_genes": 80, "pct_counts_mt": 50,
                    "pct_counts_ribo": 10, "multiple_metrics": 30, "total_removed": 210,
                },
                "n_cells": 5000, "figure_stem": qc_explore.S1_FIGURE_STEM,
            }
        }
        (art / "qc_explore.json").write_text(json.dumps(data))
        return paths

    def test_block_has_table_and_figure(self):
        with tempfile.TemporaryDirectory() as d:
            run_dir = Path(d)
            self._write_explore(run_dir)
            blocks = qc_explore.render_appendix_blocks(run_dir)
            self.assertIn("s1_rna_qc", blocks)
            block = blocks["s1_rna_qc"]
            # 4-column header
            self.assertIn("| parameter | threshold | cells removed | note |", block)
            # non-exclusive summary rows present
            self.assertIn("multiple_metrics", block)
            self.assertIn("total_removed", block)
            # marginal count surfaced
            self.assertIn("100", block)
            # relative figure reference into deliverables/figures
            self.assertIn(f"../../figures/{qc_explore.S1_FIGURE_STEM}.png", block)

    def test_render_plan_appendix_replaces_bullets(self):
        with tempfile.TemporaryDirectory() as d:
            run_dir = Path(d)
            self._write_explore(run_dir)
            blocks = qc_explore.render_appendix_blocks(run_dir)
            plan = {
                "workflow_branch": "rna_only",
                "stages": {
                    "s1_rna_qc": {"parameters": {
                        "k_mad": {"value": 5.0, "rationale": "should be hidden"}}},
                },
                "warnings": [],
            }
            md = render_plan_appendix(plan, blocks)
            self.assertIn("| parameter | threshold | cells removed | note |", md)
            # the parameter bullet for k_mad is replaced by the table
            self.assertNotIn("should be hidden", md)


class InMemoryAndRederiveTests(unittest.TestCase):
    def _rna_metrics_df(self, n: int = 300) -> pd.DataFrame:
        rng = np.random.default_rng(1)
        return pd.DataFrame({
            "total_counts": rng.integers(50, 8000, n).astype(float),
            "n_genes_by_counts": rng.integers(20, 3000, n).astype(float),
            "pct_counts_mt": rng.uniform(0, 25, n),
            "pct_counts_ribo": rng.uniform(0, 60, n),
        })

    def test_explore_rna_uses_in_memory_adata_without_h5ad(self):
        """The merged S0 job passes its loaded matrix; qc_explore must not read
        rna_ingest.h5ad and must persist the per-cell metrics parquet."""
        import anndata as ad
        import scipy.sparse as sp

        with tempfile.TemporaryDirectory() as d:
            run_dir = Path(d)
            paths = RunPaths(run_dir)
            art = paths.stage_dir("qc_explore")
            art.mkdir(parents=True, exist_ok=True)
            paths.deliv_figures.mkdir(parents=True, exist_ok=True)

            # No rna_ingest.h5ad on disk — the in-memory path must not need it.
            self.assertFalse(
                (run_dir / "internal" / "artifacts" / "s0_ingest" / "rna_ingest.h5ad").exists()
            )

            rng = np.random.default_rng(0)
            n_cells, n_genes = 200, 40
            X = sp.csr_matrix(rng.integers(0, 80, size=(n_cells, n_genes)).astype(float))
            a = ad.AnnData(X=X)
            a.var_names = [f"GENE{i}" for i in range(n_genes - 2)] + ["MT-CO1", "RPS3"]
            a.obs_names = [f"cell{i}" for i in range(n_cells)]

            plan = {"stages": {"s1_rna_qc": {"parameters": {}}}}
            res = qc_explore._explore_rna(run_dir, plan, paths.deliv_figures, art, adata=a)

            self.assertIsNotNone(res)
            self.assertEqual(res["n_cells"], n_cells)
            self.assertEqual(res["metrics_parquet"], qc_explore.RNA_METRICS_PARQUET)
            parquet = art / qc_explore.RNA_METRICS_PARQUET
            self.assertTrue(parquet.exists())
            df = pd.read_parquet(parquet)
            self.assertEqual(list(df.columns), qc_explore.RNA_METRIC_COLS)
            self.assertEqual(len(df), n_cells)

    def test_rederive_from_metrics_rebuilds_json_without_reload(self):
        """rederive_from_metrics must rebuild qc_explore.json purely from the parquet
        (no h5ad / fragment reload), matching the direct compute."""
        with tempfile.TemporaryDirectory() as d:
            run_dir = Path(d)
            paths = RunPaths(run_dir)
            art = paths.stage_dir("qc_explore")
            art.mkdir(parents=True, exist_ok=True)
            paths.deliv_figures.mkdir(parents=True, exist_ok=True)

            plan = {"stages": {"s1_rna_qc": {"parameters": {}}}}
            plan_path = paths.artifact("p2_plan", "preprocessing_plan.json")
            plan_path.parent.mkdir(parents=True, exist_ok=True)
            plan_path.write_text(json.dumps(plan))

            df = self._rna_metrics_df()
            _io.write_parquet_safe(df, art / qc_explore.RNA_METRICS_PARQUET)

            out_path = qc_explore.rederive_from_metrics(run_dir)
            data = json.loads(out_path.read_text())
            self.assertIn("s1_rna_qc", data)

            expected = qc_explore._rna_qc_from_metrics(df, {}, paths.deliv_figures)
            self.assertEqual(data["s1_rna_qc"]["cells_removed"], expected["cells_removed"])
            self.assertEqual(data["s1_rna_qc"]["n_cells"], len(df))

    def test_run_prefers_rederive_when_parquet_exists(self):
        """Standalone run() with no in-memory matrix and an existing parquet takes
        the cheap re-derive path — no rna_ingest.h5ad required."""
        with tempfile.TemporaryDirectory() as d:
            run_dir = Path(d)
            paths = RunPaths(run_dir)
            art = paths.stage_dir("qc_explore")
            art.mkdir(parents=True, exist_ok=True)
            paths.deliv_figures.mkdir(parents=True, exist_ok=True)

            plan = {"stages": {"s1_rna_qc": {"parameters": {}}}}
            plan_path = paths.artifact("p2_plan", "preprocessing_plan.json")
            plan_path.parent.mkdir(parents=True, exist_ok=True)
            plan_path.write_text(json.dumps(plan))
            _io.write_parquet_safe(self._rna_metrics_df(), art / qc_explore.RNA_METRICS_PARQUET)

            out_path = qc_explore.run(run_dir)
            data = json.loads(out_path.read_text())
            self.assertIn("s1_rna_qc", data)


if __name__ == "__main__":
    unittest.main()
