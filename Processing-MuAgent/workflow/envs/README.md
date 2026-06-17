# workflow/envs/

Environment definitions for MuAgene. **Provisioning is owned by Execution-MuAgent**
(`init-machine` / `provision-env` / `validate-env`); these files are the *definitions*
it provisions from. `manifest.yaml` is the single source of the per-device provider +
paths (read by both agents); `site.config`'s `environments:` section is generated from it.

## CPU env — `processing.yaml` (canonical `muagene` env)

Source-of-truth for the CPU env (provider `lock`: a conda-lock `processing.linux-64.lock`
generated from this YAML gives reproducible, solve-free installs). Includes:

- Python stack (scanpy, muon, snapatac2, scrublet, snakemake, …)
- **ATAC fragment prep:** `htslib` (`bgzip` + `tabix`) — must live in the env so cluster
  child jobs (which activate `$PMA_CONDA_ENV`) have them on PATH.
- **S1a ambient correction:** `r-base`, `bioconductor-celda` (DecontX), `r-soupx` (SoupX).

## GPU env — `muagene-gpu.def` (container, **pull-only**)

The GPU env is a **container**, not a conda env: the full RAPIDS + single-cell conda
solve OOMs / hits channel errors, so we bootstrap from NVIDIA's prebuilt RAPIDS base
image and add the single-cell libs with pip. The `%post` package list is the single
manifest (there is deliberately no parallel GPU YAML). This `.def` is built + pushed to
a registry **once, centrally** by a maintainer (`scripts/build_and_push_gpu_image.sh`);
every target machine then **pulls** that pinned image (`gpu_image_uri`) — no machine
builds the container locally, so there is no `--fakeroot`/subuid requirement. GPU-capable
stages run inside it via `singularity exec --nv`; `rapids-singlecell` supplies the GPU
drop-ins (`pp.scrublet`, `pp.pca`, `pp.neighbors`, `tl.leiden`, `tl.umap`, `pp.harmony_integrate`).

**Container bind contract** (`workflow/profiles/*/{slurm,pbs}-submit.sh`): the wrapper binds
**both** the resolved repo root (`PMA_REPO_ROOT`, for `launch_runner.sh` + the `executor`
package on `PYTHONPATH`) **and** the resolved run directory (`PMA_RUN_DIR`, for
`internal/artifacts/…` I/O) — the run data may sit under a nested mount that singularity's
default `$HOME`/`$PWD` auto-mount does not cover, so both are bound explicitly. An optional
extra bind (`PMA_GPU_BIND`, from `site.config` `common.scratch` / `configure-execution
--scratch`) appends after, for paths a stage writes outside the run dir.

## `*.imports.txt`

Module lists `validate-env` imports to confirm an env is usable before submitting
(`muagene.imports.txt` for CPU, `muagene-gpu.imports.txt` for GPU — the latter checks
`rapids_singlecell`/`cupy` so a CPU-only env fails loud rather than silently degrading).

## Regenerating the CPU lock

`processing.yaml` is the human source of truth; `processing.linux-64.lock` is what
actually installs. After editing the YAML, regenerate and commit the lock — otherwise
`validate-env`/`submit` fail loud with `lock_stale_vs_yaml`:

```bash
Processing-MuAgent regenerate-locks          # runs conda-lock; needs: pip install '.[dev]'
```

## Provisioning

```bash
# Execution-MuAgent owns this. Fresh machine — one bootstrap command (creates the CPU env
# from the lock, installs both packages, pulls the GPU image):
Execution-MuAgent init-machine --processing-repo <Processing-MuAgent> --device both \
  --gpu-image-uri docker://<registry>/muagene-gpu:<tag>
# Per-run (optional once bootstrapped); `submit` also auto-provisions a missing/stale env:
Execution-MuAgent provision-env --site-config <site.config> --repo-root <Processing-MuAgent> --device both
```

Whether S1a runs correction vs pass-through is set in the preprocessing plan
(`s1a_ambient.method`: default `auto`). See `executor/plan_assembler.py` and `plan_review.md`.
