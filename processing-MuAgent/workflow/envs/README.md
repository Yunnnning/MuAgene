# workflow/envs/

Reserved for per-rule Snakemake conda environment YAMLs (`scanpy_muon.yaml`,
`snapatac2.yaml`, `context.yaml`). The current MVP runs against a single shared
env (see the top-level README). Populate this directory with per-rule envs when
enabling `snakemake --use-conda` for production runs.
