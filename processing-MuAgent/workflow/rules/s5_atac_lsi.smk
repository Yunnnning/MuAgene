rule s5_atac_lsi_propose:
    input:
        # s5 consumes ATAC post-doublet (the s3 output for the ATAC side).
        atac_h5 = str(INTERNAL / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"),
    output:
        proposal = str(INTERNAL / "proposals" / "s5_atac_lsi.yaml"),
        awaiting = str(INTERNAL / "proposals" / "s5_atac_lsi.awaiting_approval"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        from executor import approval
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s5_atac_lsi",
            "action": "SnapATAC2 tile matrix + TF-IDF + spectral embedding; export peak matrix when possible, else verified tile-matrix fallback for integration",
        }))
        approval.mark_awaiting(params.run_dir, "s5_atac_lsi")


rule s5_atac_lsi_execute:
    input:
        proposal         = str(INTERNAL / "proposals" / "s5_atac_lsi.yaml"),
        approved         = str(INTERNAL / "checkpoints" / "s5_atac_lsi.approved"),
        plan             = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done = str(INTERNAL / "checkpoints" / "plan_review.approved"),
    output:
        summary = str(INTERNAL / "artifacts" / "s5_atac_lsi" / "lsi_summary.json"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import json
        from pathlib import Path
        from executor.stages import s5_atac_lsi
        plan = json.loads(Path(input.plan).read_text())
        s5_atac_lsi.run(params.run_dir, plan)
