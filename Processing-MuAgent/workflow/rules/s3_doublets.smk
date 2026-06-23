def _s3_inputs(wildcards):
    from executor import provenance
    branch = provenance.current_branch(str(INTERNAL / "parameters.yaml"))
    paths: dict = {}
    # Depend on the durable qc_summary.json marker, NOT the rna_qc.h5ad/atac_qc.h5ad
    # working files (same pattern as s2_atac_qc's _s2_propose_inputs). The marker is
    # co-written with the h5ad at stage completion and survives post_qc_review
    # cleanup, so it carries the S1/S2 -> S3 ordering edge without making the
    # deletable h5ads part of the declared DAG. S3 reads the h5ads by path.
    if branch in ("paired", "separate", "rna_only"):
        paths["rna_done"] = str(INTERNAL / "artifacts" / "s1_rna_qc" / "qc_summary.json")
    if branch in ("paired", "separate", "atac_only"):
        paths["atac_done"] = str(INTERNAL / "artifacts" / "s2_atac_qc" / "qc_summary.json")
    return paths


rule s3_doublets_propose:
    input:
        unpack(_s3_inputs)
    output:
        proposal = str(INTERNAL / "proposals" / "s3_doublets.yaml"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s3_doublets",
            "action": "Scrublet (RNA) + heuristic ATAC doublet flagging; overlap + goal-based policy recommendation",
        }))


rule s3_doublets_execute:
    input:
        unpack(_s3_inputs),
        plan             = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done = str(INTERNAL / "checkpoints" / "plan_review.approved"),
    output:
        # calls.parquet is the durable stage-done marker — written on every branch
        # and the DAG edge S3 -> qc_handoff. The post-doublet h5ads
        # (rna_post_doublet.h5ad, atac_post_doublet.h5ad) are written by the run body
        # below but are intentionally NOT declared outputs: they are transient working
        # files that qc_handoff deletes once the post-QC h5mu exists (same pattern as
        # rna_qc.h5ad vs qc_summary.json in S1/S2). Declaring them would make their
        # deletion look like a missing output and trigger a spurious S3 re-run on `all`.
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
        from executor.cluster_exit import finalize_cluster_exit
        plan = json.loads(Path(input.plan).read_text())
        branch = provenance.current_branch(str(INTERNAL / "parameters.yaml"))
        s3_doublets.run(params.run_dir, plan, workflow_branch=branch)
        finalize_cluster_exit()
