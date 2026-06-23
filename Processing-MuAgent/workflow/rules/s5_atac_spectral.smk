rule s5_atac_spectral_propose:
    input:
        # s5 rebuilds the ATAC working object from the canonical post-QC h5mu
        # (atac mod = fragments + chrom sizes), not the transient s3 h5ad.
        post_qc_h5mu = str(QC / f"post_qc_{RUN_DIR.name}.h5mu"),
    output:
        proposal = str(INTERNAL / "proposals" / "s5_atac_spectral.yaml"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s5_atac_spectral",
            "action": "SnapATAC2 tile matrix + feature selection + spectral embedding (snap.tl.spectral); export peak matrix when possible, else verified tile-matrix fallback for integration",
        }))


rule s5_atac_spectral_execute:
    input:
        plan               = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done   = str(INTERNAL / "checkpoints" / "plan_review.approved"),
        qc_review_done     = str(INTERNAL / "checkpoints" / "post_qc_review.approved"),
        post_qc_h5mu       = str(QC / f"post_qc_{RUN_DIR.name}.h5mu"),
    output:
        summary = str(INTERNAL / "artifacts" / "s5_atac_spectral" / "spectral_summary.json"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s5_atac_spectral"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s5_atac_spectral", attempt),
        runtime=RUNTIME["s5_atac_spectral"],
    run:
        import json
        from pathlib import Path
        from executor.stages import s5_atac_spectral
        from executor.cluster_exit import finalize_cluster_exit
        plan = json.loads(Path(input.plan).read_text())
        s5_atac_spectral.run(params.run_dir, plan)
        finalize_cluster_exit()
