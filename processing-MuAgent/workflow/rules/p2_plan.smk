rule p2_plan_propose:
    input:
        context = str(INTERNAL / "artifacts" / "p1_context" / "context_extraction.json"),
        ingest  = str(INTERNAL / "artifacts" / "s0_ingest" / "validation_report.json"),
    output:
        proposal = str(INTERNAL / "proposals" / "p2_plan.yaml"),
        awaiting = str(INTERNAL / "proposals" / "p2_plan.awaiting_approval"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import json
        import yaml
        from pathlib import Path
        from executor import plan_assembler as pa, approval, provenance
        ctx = json.loads(Path(input.context).read_text())
        # Sample type from context extraction
        sample_type = (ctx.get("fields", {}).get("sample_type") or {}).get("value", "unknown")
        # Workflow branch from S0
        branch = provenance.get_value(str(INTERNAL / "parameters.yaml"),
                                      "plan.workflow_branch", "paired")
        study_goal = config.get("study_goal")
        plan = pa.assemble_plan(params.run_dir,
                                workflow_branch=branch,
                                sample_type=sample_type,
                                study_goal=study_goal)
        _, plan_hash = pa.write_plan(params.run_dir, plan)
        # Merged plan_review.md (summary + appendix) is written at plan_review_propose.
        provenance.set_param(
            str(INTERNAL / "parameters.yaml"), "plan.plan_hash", plan_hash,
            source="derived", confidence="high",
            rationale="sha256 of preprocessing_plan.json",
            method={"name": "sha256_bytes", "code_ref": "executor/hashing.py::sha256_bytes"},
        )
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "p2_plan",
            "branch": branch,
            "sample_type": sample_type,
            "plan_hash": plan_hash,
        }))
        approval.mark_awaiting(params.run_dir, "p2_plan")


rule p2_plan_execute:
    input:
        proposal = str(INTERNAL / "proposals" / "p2_plan.yaml"),
        approved = str(INTERNAL / "checkpoints" / "p2_plan.approved"),
    output:
        plan     = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
    run:
        # Plan was written at propose time; re-assemble deterministically so execute has
        # a real output (plan is pure function of context + ingest + config).
        import json
        from pathlib import Path
        from executor import plan_assembler as pa, provenance
        ctx_path = Path(input.proposal).parent.parent / "artifacts" / "p1_context" / "context_extraction.json"
        ctx = json.loads(ctx_path.read_text()) if ctx_path.exists() else {"fields": {}}
        sample_type = (ctx.get("fields", {}).get("sample_type") or {}).get("value", "unknown")
        branch = provenance.get_value(str(INTERNAL / "parameters.yaml"),
                                      "plan.workflow_branch", "paired")
        study_goal = config.get("study_goal")
        plan = pa.assemble_plan(str(RUN_DIR), workflow_branch=branch,
                                sample_type=sample_type, study_goal=study_goal)
        out_path, plan_hash = pa.write_plan(str(RUN_DIR), plan)
        provenance.set_param(
            str(INTERNAL / "parameters.yaml"), "plan.plan_hash", plan_hash,
            source="derived", confidence="high",
            rationale="sha256 of preprocessing_plan.json",
            method={"name": "sha256_bytes", "code_ref": "executor/hashing.py::sha256_bytes"},
        )
