rule post_qc_review_propose:
    """QC review user checkpoint (#2): figures + qc_review_<run>.md after S3.

    Generates:
      - deliverables/checkpoint/qc_review/figures/post_qc_review_cell_counts.{png,pdf}
      - deliverables/checkpoint/qc_review/figures/post_qc_review_doublet_rna.{png,pdf}
      - deliverables/checkpoint/qc_review/figures/post_qc_review_doublet_atac.{png,pdf}
      - deliverables/checkpoint/qc_review/qc_review_<run>.md
      - deliverables/checkpoint/qc_review/qc_summary_<run>.html

    S1/S2 QC figures are already in deliverables/checkpoint/qc_review/figures/.
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
        summary  = str(CHECKPOINT / "qc_review" / f"qc_review_{RUN_DIR.name}.md"),
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
                "QC review checkpoint — inspect deliverables/checkpoint/qc_review/figures/ "
                f"and qc_review_{RUN_DIR.name}.md (includes paired S3 union doublet policy when applicable). "
                "Revise S1/S2 thresholds if needed, "
                "then run: Processing-MuAgent approve post_qc_review --config $CFG"
            ),
            **result,
        }))
        approval.mark_awaiting(params.run_dir, "post_qc_review")
