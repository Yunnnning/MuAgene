"""Tests for S1a marker gene resolution and plotting."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from click.testing import CliRunner

from executor import cli
from executor import provenance as _prov
from executor.stages import s1a_ambient as s1a


def _make_adata(n_obs: int = 10, genes: list[str] | None = None) -> ad.AnnData:
    genes = genes or ["Kit", "Stra8", "Sycp3", "Meiob", "Acrv1"]
    rng = np.random.default_rng(0)
    raw = sp.csr_matrix(rng.integers(0, 20, size=(n_obs, len(genes)), dtype=np.int32))
    corr = sp.csr_matrix(rng.integers(0, 15, size=(n_obs, len(genes)), dtype=np.int32))
    obs_names = [f"cell_{i}" for i in range(n_obs)]
    return ad.AnnData(
        X=corr,
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=genes),
        layers={"counts_raw": raw, "counts": corr},
    )


def _write_run_layout(
    tmp: Path,
    *,
    genes: list[str] | None = None,
    yaml_genes: list[str] | None = None,
) -> Path:
    art = tmp / "internal" / "artifacts" / "s1a_ambient"
    art.mkdir(parents=True, exist_ok=True)
    params_path = tmp / "internal" / "parameters.yaml"
    params_path.parent.mkdir(parents=True, exist_ok=True)
    if not params_path.exists():
        params_path.write_text("{}\n")

    a = _make_adata(genes=genes)
    h5ad_path = art / "rna_decontaminated.h5ad"
    a.write_h5ad(h5ad_path)

    coords = np.column_stack([
        np.linspace(0, 1, a.n_obs),
        np.linspace(1, 0, a.n_obs),
    ])
    pd.DataFrame({
        "obs_name": list(a.obs_names),
        "tsne_x": coords[:, 0].astype("float32"),
        "tsne_y": coords[:, 1].astype("float32"),
    }).to_parquet(art / "tsne_coords_cache.parquet", index=False)

    raw = np.asarray(a.layers["counts_raw"].sum(axis=1)).ravel()
    corr = np.asarray(a.layers["counts"].sum(axis=1)).ravel()
    pd.DataFrame({
        "obs_name": list(a.obs_names),
        "total_raw": raw.astype("float32"),
        "total_corrected": corr.astype("float32"),
    }).to_parquet(art / "cell_totals.parquet", index=False)

    if yaml_genes is not None:
        _prov.set_param(
            params_path, "s1a_ambient.marker_genes", yaml_genes,
            source="user", confidence="high",
            rationale="test fixture",
        )

    cfg = tmp / "run.yaml"
    cfg.write_text(yaml.safe_dump({"run_dir": str(tmp)}))
    return tmp


class ResolveMarkerGenesTests(unittest.TestCase):
    def test_yaml_overrides_plan(self):
        with tempfile.TemporaryDirectory() as tmp:
            params_path = Path(tmp) / "parameters.yaml"
            params_path.write_text("{}\n")
            _prov.set_param(
                params_path, "s1a_ambient.marker_genes", ["Kit", "Stra8"],
                source="user", confidence="high", rationale="test",
            )
            plan = {"marker_genes": {"value": ["Other"]}}
            self.assertEqual(
                s1a.resolve_marker_genes(params_path, plan),
                ["Kit", "Stra8"],
            )

    def test_plan_fallback_when_yaml_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            params_path = Path(tmp) / "parameters.yaml"
            params_path.write_text("{}\n")
            plan = {"marker_genes": {"value": ["Kit", "Meiob"]}}
            self.assertEqual(
                s1a.resolve_marker_genes(params_path, plan),
                ["Kit", "Meiob"],
            )

    def test_empty_when_both_unset(self):
        with tempfile.TemporaryDirectory() as tmp:
            params_path = Path(tmp) / "parameters.yaml"
            params_path.write_text("{}\n")
            self.assertEqual(s1a.resolve_marker_genes(params_path, {}), [])


class MatchGenesTests(unittest.TestCase):
    def test_case_insensitive_match(self):
        a = _make_adata()
        found, missing = s1a._match_genes_to_var(a, ["kit", "STRA8"])
        self.assertEqual(found, ["Kit", "Stra8"])
        self.assertEqual(missing, [])

    def test_preserves_request_order(self):
        a = _make_adata()
        found, missing = s1a._match_genes_to_var(
            a, ["Meiob", "Kit", "Sycp3"],
        )
        self.assertEqual(found, ["Meiob", "Kit", "Sycp3"])
        self.assertEqual(missing, [])


class PlotMarkerGenesTests(unittest.TestCase):
    def test_cache_hit_uses_backed_h5ad_without_tsne_recompute(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = _write_run_layout(Path(tmp))
            tsne_calls: list[str] = []

            def _fail_tsne(*args, **kwargs):
                tsne_calls.append("called")
                raise AssertionError("t-SNE should not run on cache hit")

            with mock.patch.object(s1a, "_load_or_compute_tsne", side_effect=_fail_tsne):
                with mock.patch("executor.figures.plot_marker_genes_tsne") as plot_fn:
                    result = s1a._plot_marker_genes(
                        run_dir, ["Kit", "Stra8"],
                        write_params=False, refresh_qc=False,
                    )
            self.assertEqual(tsne_calls, [])
            self.assertEqual(result["found"], ["Kit", "Stra8"])
            plot_fn.assert_called_once()
            self.assertEqual(plot_fn.call_args.kwargs["genes"], ["Kit", "Stra8"])

    def test_gene_order_passed_to_plot(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = _write_run_layout(Path(tmp))
            order = ["Acrv1", "Meiob", "Kit"]
            with mock.patch("executor.figures.plot_marker_genes_tsne") as plot_fn:
                s1a._plot_marker_genes(
                    run_dir, order, write_params=False, refresh_qc=False,
                )
            self.assertEqual(plot_fn.call_args.kwargs["genes"], order)

    def test_writes_marker_check_artifact(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = _write_run_layout(Path(tmp))
            with mock.patch("executor.figures.plot_marker_genes_tsne"):
                s1a._plot_marker_genes(
                    run_dir, ["Kit", "NotReal"],
                    write_params=False, refresh_qc=False,
                )
            art = run_dir / "internal" / "artifacts" / "s1a_ambient" / "marker_gene_check.json"
            data = json.loads(art.read_text())
            self.assertEqual(data["found"], ["Kit"])
            self.assertEqual(data["missing"], ["NotReal"])


class MarkerGeneCheckCliTests(unittest.TestCase):
    def test_empty_found_does_not_claim_figure_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = _write_run_layout(Path(tmp))
            cfg = run_dir / "run.yaml"
            with mock.patch.object(
                s1a, "run_marker_gene_check",
                return_value={"found": [], "missing": ["FakeGene"]},
            ):
                runner = CliRunner()
                result = runner.invoke(
                    cli.main,
                    ["marker-gene-check", "--config", str(cfg), "FakeGene"],
                )
            self.assertEqual(result.exit_code, 0)
            self.assertIn("No marker genes found in matrix; figure not written.", result.output)
            self.assertNotIn("Figure written", result.output)
            self.assertNotIn("QC reports refreshed", result.output)

    def test_found_refreshes_qc_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = _write_run_layout(Path(tmp))
            cfg = run_dir / "run.yaml"
            with mock.patch.object(
                s1a, "run_marker_gene_check",
                return_value={"found": ["Kit"], "missing": []},
            ) as check_fn:
                runner = CliRunner()
                result = runner.invoke(
                    cli.main,
                    ["marker-gene-check", "--config", str(cfg), "Kit"],
                )
            self.assertEqual(result.exit_code, 0)
            check_fn.assert_called_once()
            self.assertTrue(check_fn.call_args.kwargs.get("refresh_qc"))
            self.assertIn("QC reports refreshed.", result.output)

    def test_plot_only_skips_qc_refresh(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = _write_run_layout(Path(tmp))
            cfg = run_dir / "run.yaml"
            with mock.patch.object(
                s1a, "run_marker_gene_check",
                return_value={"found": ["Kit"], "missing": []},
            ) as check_fn:
                runner = CliRunner()
                result = runner.invoke(
                    cli.main,
                    [
                        "marker-gene-check", "--plot-only",
                        "--config", str(cfg), "Kit",
                    ],
                )
            self.assertEqual(result.exit_code, 0)
            self.assertFalse(check_fn.call_args.kwargs.get("refresh_qc"))
            self.assertIn("--plot-only: QC reports unchanged", result.output)


if __name__ == "__main__":
    unittest.main()
