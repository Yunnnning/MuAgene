rule s7_clustering_propose:
    output:
        proposal = str(INTERNAL / "proposals" / "s7_clustering.yaml"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s7_clustering",
            "action": "Leiden clustering at fixed resolutions (RNA=0.7, ATAC=0.5); "
                      "no sweep, no checkpoint",
        }))


rule s7_clustering_execute:
    """Cluster each modality at its fixed resolution and write final labels.
    Runs automatically after S6 — there is no resolution-review checkpoint.
    """
    input:
        plan             = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done = str(INTERNAL / "checkpoints" / "plan_review.approved"),
        # Durable S6 marker, not the deletable rna_neighbors.h5ad (read by path in-stage).
        rna_neighbors    = str(INTERNAL / "artifacts" / "s6_neighbors" / "neighbors_summary.json"),
    output:
        # clustering_summary.json is the SOLE declared output + durable stage-done
        # marker (S7 -> S8 edge). rna_clustered.h5ad and atac_leiden_labels.parquet are
        # UNTRACKED working files read by S8 by path and removed by `finish-cleanup`.
        summary = str(INTERNAL / "artifacts" / "s7_clustering" / "clustering_summary.json"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s7_clustering"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s7_clustering", attempt),
        runtime=RUNTIME["s7_clustering"],
    run:
        import json
        from pathlib import Path
        from executor.stages import s7_clustering
        from executor.cluster_exit import finalize_cluster_exit
        plan = json.loads(Path(input.plan).read_text())
        s7_clustering.run(params.run_dir, plan)
        finalize_cluster_exit()
