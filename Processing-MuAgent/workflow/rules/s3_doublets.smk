def _s3_inputs(wildcards):
    from executor import provenance
    branch = provenance.current_branch(str(INTERNAL / "parameters.yaml"))
    paths: dict = {}
    if branch in ("paired", "separate", "rna_only"):
        paths["rna_qc"] = str(INTERNAL / "artifacts" / "s1_rna_qc" / "rna_qc.h5ad")
    if branch in ("paired", "separate", "atac_only"):
        paths["atac_qc"] = str(INTERNAL / "artifacts" / "s2_atac_qc" / "atac_qc.h5ad")
    return paths


rule s3_doublets_propose:
    input:
        unpack(_s3_inputs)
    output:
        proposal = str(INTERNAL / "proposals" / "s3_doublets.yaml"),
        awaiting = str(INTERNAL / "proposals" / "s3_doublets.awaiting_approval"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        from executor import approval
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s3_doublets",
            "action": "Scrublet (RNA) + heuristic ATAC doublet flagging; overlap + goal-based policy recommendation",
        }))
        approval.mark_awaiting(params.run_dir, "s3_doublets")


rule s3_doublets_execute:
    input:
        proposal         = str(INTERNAL / "proposals" / "s3_doublets.yaml"),
        approved         = str(INTERNAL / "checkpoints" / "s3_doublets.approved"),
        plan             = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done = str(INTERNAL / "checkpoints" / "plan_review.approved"),
    output:
        # Both h5ad outputs are always declared; the stage writes an empty
        # placeholder for the missing-modality side in single-modality branches
        # so the declared DAG is satisfied regardless of branch.
        rna_post  = str(INTERNAL / "artifacts" / "s3_doublets" / "rna_post_doublet.h5ad"),
        atac_post = str(INTERNAL / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"),
        calls     = str(INTERNAL / "artifacts" / "s3_doublets" / "calls.parquet"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s3_doublets"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s3_doublets", attempt),
        runtime=RUNTIME["s3_doublets"],
    run:
        import json
        from pathlib import Path
        from executor.stages import s3_doublets
        from executor import provenance
        plan = json.loads(Path(input.plan).read_text())
        branch = provenance.current_branch(str(INTERNAL / "parameters.yaml"))
        s3_doublets.run(params.run_dir, plan, workflow_branch=branch)
