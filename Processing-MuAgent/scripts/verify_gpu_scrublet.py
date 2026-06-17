#!/usr/bin/env python
"""Phase-0 POC: verify CPU/GPU Scrublet parity before wiring GPU into the executor.

Runs CPU ``scrublet.Scrublet`` and GPU ``rapids_singlecell.pp.scrublet`` on the SAME
count matrix and checks that the doublet scores agree (rank correlation) and the calls
concordance, plus timing. Pins the exact rapids API used by
``executor/stages/s3_doublets.py:_rna_scrublet_gpu``.

Run INSIDE the GPU container on a GPU node, e.g.:

  module load singularityce/3.11.3
  srun -p gpu --gres=gpu:A5000:1 --account=vaquerizas \
    singularity exec --nv ~/.muagene/images/muagene-gpu.sif \
    python scripts/verify_gpu_scrublet.py \
      [--h5ad <run>/internal/artifacts/s1_rna_qc/rna_qc.h5ad]

With no --h5ad it synthesizes a small count matrix so the script is self-contained.
Exits non-zero if parity is below tolerance.
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import scipy.sparse as sp


def _load_counts(h5ad: str | None):
    if h5ad:
        import anndata as ad
        a = ad.read_h5ad(h5ad)
        X = a.layers["counts"] if "counts" in a.layers else a.X
        return sp.csr_matrix(X), f"{h5ad} ({a.n_obs} cells x {a.n_vars} genes)"
    # Synthetic: NB-ish counts + a doublet population (summed pairs).
    rng = np.random.default_rng(0)
    n, g = 2000, 1500
    base = rng.poisson(0.2, size=(n, g)).astype(np.float32)
    dbl = base[rng.integers(0, n, 200)] + base[rng.integers(0, n, 200)]
    X = np.vstack([base, dbl])
    return sp.csr_matrix(X), f"synthetic ({X.shape[0]} cells x {X.shape[1]} genes)"


def _cpu(counts, rate):
    import scrublet as scr
    t = time.time()
    sd = scr.Scrublet(counts, expected_doublet_rate=rate, random_state=0)
    scores, _ = sd.scrub_doublets(verbose=False)
    return np.asarray(scores, dtype=float), time.time() - t


def _gpu(counts, rate):
    import anndata as ad
    import rapids_singlecell as rsc
    t = time.time()
    a = ad.AnnData(X=counts.copy())
    mover = getattr(getattr(rsc, "get", None), "anndata_to_GPU", None) \
        or getattr(getattr(rsc, "utils", None), "anndata_to_GPU", None)
    if mover is not None:
        mover(a)
    rsc.pp.scrublet(a, expected_doublet_rate=rate, random_state=0)
    col = "doublet_score" if "doublet_score" in a.obs else "scrublet_score"
    return np.asarray(a.obs[col].to_numpy(), dtype=float), time.time() - t


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=None)
    ap.add_argument("--rate", type=float, default=0.06)
    ap.add_argument("--threshold", type=float, default=0.25)
    ap.add_argument("--min-corr", type=float, default=0.9)
    ap.add_argument("--min-concordance", type=float, default=0.9)
    args = ap.parse_args()

    counts, desc = _load_counts(args.h5ad)
    print(f"data: {desc}")

    cpu_scores, cpu_t = _cpu(counts, args.rate)
    print(f"CPU  scrublet.Scrublet:        {cpu_t:6.1f}s")
    gpu_scores, gpu_t = _gpu(counts, args.rate)
    print(f"GPU  rapids_singlecell.scrublet:{gpu_t:6.1f}s  (speedup {cpu_t / max(gpu_t, 1e-6):.1f}x)")

    from scipy.stats import spearmanr
    corr = float(spearmanr(cpu_scores, gpu_scores).correlation)
    cpu_calls = cpu_scores > args.threshold
    gpu_calls = gpu_scores > args.threshold
    concordance = float((cpu_calls == gpu_calls).mean())
    print(f"score rank-correlation: {corr:.3f} (min {args.min_corr})")
    print(f"call concordance @ {args.threshold}: {concordance:.3f} (min {args.min_concordance}); "
          f"CPU {cpu_calls.sum()} vs GPU {gpu_calls.sum()} doublets")

    ok = corr >= args.min_corr and concordance >= args.min_concordance
    print("PARITY: PASS" if ok else "PARITY: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
