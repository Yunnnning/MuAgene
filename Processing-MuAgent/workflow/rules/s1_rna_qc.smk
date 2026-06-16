rule s1_rna_qc_propose:
    input:
        plan              = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        rna               = str(INTERNAL / "artifacts" / "s1a_ambient" / "rna_decontaminated.h5ad"),
        plan_review_done  = str(INTERNAL / "checkpoints" / "plan_review.approved"),
        plan_review_md    = str(PLAN / f"plan_review_{RUN_DIR.name}.md"),
    output:
        proposal = str(INTERNAL / "proposals" / "s1_rna_qc.yaml"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s1_rna_qc",
            "action": "MAD thresholds on total_counts, n_genes, pct_counts_mt",
        }))


rule s1_rna_qc_execute:
    input:
        plan              = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done  = str(INTERNAL / "checkpoints" / "plan_review.approved"),
        rna_decontaminated = str(INTERNAL / "artifacts" / "s1a_ambient" / "rna_decontaminated.h5ad"),
    output:
        # qc_summary.json is the SOLE declared output and the durable stage-done
        # marker (status, reports, and the s3 dependency edge all key off it).
        # rna_qc.h5ad is written by the stage as an UNTRACKED working file: it is
        # consumed only by s3_doublets (read by path) and removed by
        # _cleanup_qc_intermediates at post_qc_review approval. Keeping it out of
        # the declared DAG means deleting it never triggers a "Missing output
        # files" re-run of S1/S3.
        qc_summary = str(INTERNAL / "artifacts" / "s1_rna_qc" / "qc_summary.json"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s1_rna_qc"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s1_rna_qc", attempt),
        runtime=RUNTIME["s1_rna_qc"],
    run:
        import json
        from pathlib import Path
        from executor.stages import s1_rna_qc
        from executor.cluster_exit import finalize_cluster_exit
        plan = json.loads(Path(input.plan).read_text())
        s1_rna_qc.run(params.run_dir, plan)
        finalize_cluster_exit()
