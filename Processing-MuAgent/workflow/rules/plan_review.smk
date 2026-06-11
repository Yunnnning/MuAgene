rule plan_review_propose:
    """Render the merged plan review summary (concise + appendix). Writes the awaiting_approval sentinel.
    User must approve before any preprocessing stage (S1..S8) runs.
    """
    input:
        plan        = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        ingest      = str(INTERNAL / "artifacts" / "s0_ingest" / "validation_report.json"),
        ctx         = str(INTERNAL / "artifacts" / "p1_context" / "context_extraction.json"),
        qc_explore  = str(INTERNAL / "artifacts" / "qc_explore" / "qc_explore.json"),
    output:
        # plan_review.md is a pre-run deliverable (reviewed before plan approval).
        summary      = str(PLAN / "summary" / "plan_review.md"),
        summary_html = str(PLAN / "summary" / "plan_summary.html"),
        awaiting     = str(INTERNAL / "proposals" / "plan_review.awaiting_approval"),
    params:
        run_dir = str(RUN_DIR),
    run:
        from executor import plan_review as pr, approval
        out = pr.write_summary(params.run_dir)
        pr.write_plan_summary_html(params.run_dir)
        approval.mark_awaiting(params.run_dir, "plan_review")
