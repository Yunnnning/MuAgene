rule post_qc_review_propose:
    """QC review user checkpoint (#2): figures + qc_review.md after S3.

    Generates:
      - deliverables/checkpoint/qc_review/post_qc_review_cell_counts.{png,pdf}
      - deliverables/checkpoint/qc_review/post_qc_review_doublet_rna.{png,pdf}
      - deliverables/checkpoint/qc_review/post_qc_review_doublet_atac.{png,pdf}
      - deliverables/checkpoint/qc_review/qc_review.md

    S1/S2 QC figures are already in deliverables/checkpoint/qc_review/.
    On paired multiome, the summary documents the S3 union doublet policy for
    confirmation at this checkpoint.
    """
    input:
        rna_post  = str(INTERNAL / "artifacts" / "s3_doublets" / "rna_post_doublet.h5ad"),
        atac_post = str(INTERNAL / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"),
        calls     = str(INTERNAL / "artifacts" / "s3_doublets" / "calls.parquet"),
    output:
        proposal = str(INTERNAL / "proposals" / "post_qc_review.yaml"),
        awaiting = str(INTERNAL / "proposals" / "post_qc_review.awaiting_approval"),
        summary  = str(CHECKPOINT / "qc_review" / "qc_review.md"),
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
                "QC review checkpoint — inspect deliverables/checkpoint/qc_review/ "
                "and qc_review.md (includes paired S3 union doublet policy when applicable). "
                "Revise S1/S2 thresholds if needed, "
                "then run: Processing-MuAgent approve post_qc_review --config $CFG"
            ),
            **result,
        }))
        approval.mark_awaiting(params.run_dir, "post_qc_review")
