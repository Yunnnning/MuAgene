def _s6_inputs(wildcards):
    from executor import provenance
    branch = provenance.current_branch(str(INTERNAL / "parameters.yaml"))
    paths: dict = {}
    # Depend on the durable per-stage markers (norm_summary.json / spectral_summary.json),
    # NOT the rna_norm.h5ad / atac_spectral.h5ad working files that `finish-cleanup`
    # deletes. The markers carry the ordering edge and survive cleanup; the stages read
    # the h5ads by path.
    if branch in ("paired", "separate", "rna_only"):
        paths["rna_norm"] = str(INTERNAL / "artifacts" / "s4_rna_norm" / "norm_summary.json")
    if branch in ("paired", "separate", "atac_only"):
        paths["atac_spectral"] = str(INTERNAL / "artifacts" / "s5_atac_spectral" / "spectral_summary.json")
    return paths


rule s6_neighbors_propose:
    input:
        unpack(_s6_inputs)
    output:
        proposal = str(INTERNAL / "proposals" / "s6_neighbors.yaml"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s6_neighbors",
            "action": "PCA (RNA) + neighbor graph; ATAC KNN on spectral embedding (X_spectral)",
        }))


rule s6_neighbors_execute:
    input:
        unpack(_s6_inputs),
        plan             = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done = str(INTERNAL / "checkpoints" / "plan_review.approved"),
    output:
        # neighbors_summary.json is the SOLE declared output + durable stage-done
        # marker (S6 -> S7 edge). rna_neighbors.h5ad is an UNTRACKED working file read
        # by S7 by path and removed by `finish-cleanup`.
        summary = str(INTERNAL / "artifacts" / "s6_neighbors" / "neighbors_summary.json"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s6_neighbors"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s6_neighbors", attempt),
        runtime=RUNTIME["s6_neighbors"],
    run:
        import json
        from pathlib import Path
        from executor.stages import s6_neighbors
        from executor.cluster_exit import finalize_cluster_exit
        plan = json.loads(Path(input.plan).read_text())
        s6_neighbors.run(params.run_dir, plan)
        finalize_cluster_exit()
