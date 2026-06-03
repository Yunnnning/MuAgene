rule s4_rna_norm_propose:
    input:
        rna_post = str(INTERNAL / "artifacts" / "s3_doublets" / "rna_post_doublet.h5ad"),
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
        rna_post           = str(INTERNAL / "artifacts" / "s3_doublets" / "rna_post_doublet.h5ad"),
    output:
        h5ad = str(INTERNAL / "artifacts" / "s4_rna_norm" / "rna_norm.h5ad"),
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
