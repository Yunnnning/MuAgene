rule post_qc_review_propose:
    """Generate QC figures + early QC summary after S3; pause for user review
    before dimensionality reduction (S4/S5) is allowed to proceed.

    Generates:
      - deliverables/post_run/figures/post_qc_review_cell_counts.{png,pdf}
      - deliverables/post_run/figures/post_qc_review_doublet_rna.{png,pdf}
      - deliverables/post_run/figures/post_qc_review_doublet_atac.{png,pdf}
      - deliverables/post_run/summary/qc_summary_pre_dimred.md

    Note: RNA QC violin and ATAC fragment-size figures written by S1/S2 are
    also available in deliverables/post_run/figures/ at this point.
    """
    input:
        rna_post  = str(INTERNAL / "artifacts" / "s3_doublets" / "rna_post_doublet.h5ad"),
        atac_post = str(INTERNAL / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"),
        calls     = str(INTERNAL / "artifacts" / "s3_doublets" / "calls.parquet"),
    output:
        proposal = str(INTERNAL / "proposals" / "post_qc_review.yaml"),
        awaiting = str(INTERNAL / "proposals" / "post_qc_review.awaiting_approval"),
        summary  = str(POST_RUN / "summary" / "qc_summary_pre_dimred.md"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        from executor import approval
        from executor.stages import post_qc_review
        result = post_qc_review.propose(params.run_dir)
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "post_qc_review",
            "action": (
                "QC figures and summary written — review before approving. "
                "See deliverables/post_run/summary/qc_summary_pre_dimred.md "
                "and deliverables/post_run/figures/ then run: "
                "executor approve post_qc_review --config $CFG"
            ),
            **result,
        }))
        approval.mark_awaiting(params.run_dir, "post_qc_review")
