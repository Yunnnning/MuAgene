"""Tests for s2_atac_qc._acquire_atac reuse-vs-fresh-import branch selection.

The heavy SnapATAC2 fragment import is mocked: these tests assert which branch
runs (reuse of qc_explore's artifacts vs a fresh import) and that the reuse path
returns the per-cell metric arrays from the explore parquet without re-importing.
"""
from __future__ import annotations

import json
import sys
import tempfile
import types
import unittest
from contextlib import contextmanager
from pathlib import Path

import pandas as pd

from executor.stages import s2_atac_qc


class _FakeAdata:
    def __init__(self, n_obs: int) -> None:
        self.n_obs = n_obs
        self.closed = False

    def close(self) -> None:
        self.closed = True


def _fake_snap(n_obs: int) -> types.ModuleType:
    mod = types.ModuleType("snapatac2")
    mod.genome = types.SimpleNamespace(mm10=object())
    mod.read = lambda path: _FakeAdata(n_obs)
    return mod


@contextmanager
def _patched(snap_mod, fresh_result):
    prev_snap = sys.modules.get("snapatac2")
    prev_fresh = s2_atac_qc._import_atac_fresh
    calls = {"fresh": 0}

    def _fresh(*_a, **_k):
        calls["fresh"] += 1
        return fresh_result

    sys.modules["snapatac2"] = snap_mod
    s2_atac_qc._import_atac_fresh = _fresh
    try:
        yield calls
    finally:
        s2_atac_qc._import_atac_fresh = prev_fresh
        if prev_snap is None:
            sys.modules.pop("snapatac2", None)
        else:
            sys.modules["snapatac2"] = prev_snap


def _seed_explore(run_dir: Path, *, n_cells: int) -> Path:
    explore = run_dir / "internal" / "artifacts" / "qc_explore"
    explore.mkdir(parents=True)
    m = pd.DataFrame({
        "n_fragment": [float(i + 1) for i in range(n_cells)],
        "tsse": [float(i + 10) for i in range(n_cells)],
        "nucleosome_signal": [0.1 * (i + 1) for i in range(n_cells)],
    })
    m.to_parquet(explore / "atac_qc_metrics.parquet")
    h5 = explore / "atac_snap_explore.h5ad"
    h5.write_bytes(b"stub")
    (explore / "atac_explore_meta.json").write_text(json.dumps({
        "atac_snap_h5ad": str(h5),
        "metrics_parquet": "atac_qc_metrics.parquet",
        "genome": "mm10",
        "add_chr_prefix": True,
        "n_cells": n_cells,
    }))
    return explore


class AcquireAtacTests(unittest.TestCase):
    def test_reuse_path_skips_fresh_import(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            _seed_explore(run_dir, n_cells=3)
            art = run_dir / "internal" / "artifacts" / "s2_atac_qc"
            art.mkdir(parents=True)
            params = run_dir / "internal" / "parameters.yaml"

            with _patched(_fake_snap(3), ("FRESH",)) as calls:
                adata, n_frag, tss, ns, genome_ref, add_chr = s2_atac_qc._acquire_atac(
                    run_dir, art, params,
                )

            self.assertEqual(calls["fresh"], 0)
            self.assertEqual(adata.n_obs, 3)
            self.assertEqual(list(n_frag), [1.0, 2.0, 3.0])
            self.assertEqual(list(tss), [10.0, 11.0, 12.0])
            self.assertEqual([round(x, 1) for x in ns], [0.1, 0.2, 0.3])
            self.assertTrue(add_chr)
            self.assertIsNotNone(genome_ref)

    def test_fallback_when_meta_absent(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            art = run_dir / "internal" / "artifacts" / "s2_atac_qc"
            art.mkdir(parents=True)
            params = run_dir / "internal" / "parameters.yaml"

            with _patched(_fake_snap(3), ("FRESH",)) as calls:
                out = s2_atac_qc._acquire_atac(run_dir, art, params)

            self.assertEqual(calls["fresh"], 1)
            self.assertEqual(out, ("FRESH",))

    def test_fallback_on_size_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            _seed_explore(run_dir, n_cells=3)
            art = run_dir / "internal" / "artifacts" / "s2_atac_qc"
            art.mkdir(parents=True)
            params = run_dir / "internal" / "parameters.yaml"

            # adata reports 5 obs but the parquet only has 3 rows -> mismatch.
            with _patched(_fake_snap(5), ("FRESH",)) as calls:
                out = s2_atac_qc._acquire_atac(run_dir, art, params)

            self.assertEqual(calls["fresh"], 1)
            self.assertEqual(out, ("FRESH",))


if __name__ == "__main__":
    unittest.main()
