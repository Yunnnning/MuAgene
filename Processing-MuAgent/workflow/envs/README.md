# workflow/envs/

Conda environment specs for Processing-MuAgent.

## `processing.yaml` (canonical `grn` env)

Single shared env for local and HPC runs. Includes:

- Python stack (scanpy, muon, snapatac2, snakemake, …)
- **S1a ambient correction:** `r-base`, `bioconductor-celda` (DecontX), `r-soupx` (SoupX)

Recreate on a fresh site:

```bash
micromamba env create -n grn -f workflow/envs/processing.yaml
micromamba activate grn
pip install -e .   # from Processing-MuAgent root
```

Whether S1a runs correction vs pass-through is set in the preprocessing plan (`s1a_ambient.method`: default `auto`, confirm at plan review from `study_goal` and inputs; override with `s1a_ambient_method` in `run.yaml`). See `executor/plan_assembler.py` and `plan_review.md`.

Per-rule Snakemake conda env YAMLs may be added later when enabling `snakemake --use-conda` for production.
