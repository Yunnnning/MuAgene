rule s4_rna_norm_propose:
    input:
        # Reads RNA from the canonical post-QC h5mu, not the transient s3 h5ad.
        post_qc_h5mu = str(QC / f"post_qc_{RUN_DIR.name}.h5mu"),
    output:
        proposal = str(INTERNAL / "proposals" / "s4_rna_norm.yaml"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s4_rna_norm",
            "action": "log-normalize (target_sum=1e4) + HVG seurat_v3 on counts layer",
        }))


rule s4_rna_norm_execute:
    input:
        plan               = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done   = str(INTERNAL / "checkpoints" / "plan_review.approved"),
        qc_review_done     = str(INTERNAL / "checkpoints" / "post_qc_review.approved"),
        post_qc_h5mu       = str(QC / f"post_qc_{RUN_DIR.name}.h5mu"),
    output:
        # norm_summary.json is the SOLE declared output and the durable stage-done
        # marker (status + the S4 -> S6 dependency edge key off it). rna_norm.h5ad is
        # written as an UNTRACKED working file: read by S6 by path, removed by
        # `finish-cleanup` once the processed deliverable exists. Keeping it out of
        # the declared DAG means deleting it never triggers a "Missing output" re-run.
        summary = str(INTERNAL / "artifacts" / "s4_rna_norm" / "norm_summary.json"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s4_rna_norm"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s4_rna_norm", attempt),
        runtime=RUNTIME["s4_rna_norm"],
    run:
        import json
        from pathlib import Path
        from executor.stages import s4_rna_norm
        from executor.cluster_exit import finalize_cluster_exit
        plan = json.loads(Path(input.plan).read_text())
        s4_rna_norm.run(params.run_dir, plan)
        finalize_cluster_exit()
