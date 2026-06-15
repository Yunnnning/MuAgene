rule plan_review_propose:
    """Arm the plan_review gate and write per-stage job specs.

    Deliverables (plan_review_<run>.md and plan_summary_<run>.html) are written
    by the agent via ``executor plan-review --intro "..."``, not here.
    User must approve before any preprocessing stage (S1..S8) runs.
    """
    input:
        plan        = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        ingest      = str(INTERNAL / "artifacts" / "s0_ingest" / "validation_report.json"),
        ctx         = str(INTERNAL / "artifacts" / "p1_context" / "context_extraction.json"),
        qc_explore  = str(INTERNAL / "artifacts" / "qc_explore" / "qc_explore.json"),
    output:
        awaiting    = str(INTERNAL / "proposals" / "plan_review.awaiting_approval"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import json as _json
        from executor import approval, specs as _specs
        branch = _json.loads(open(input.plan).read()).get("workflow_branch", "paired")
        _specs.write_stage_specs(params.run_dir, branch)
        approval.mark_awaiting(params.run_dir, "plan_review")
