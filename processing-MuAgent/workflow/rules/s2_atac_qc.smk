def _s2_propose_inputs(wildcards):
    from executor import provenance
    branch = provenance.current_branch(str(INTERNAL / "parameters.yaml"))
    paths: dict = {
        "plan": str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
    }
    # Order-only dependency on S1 in branches where S1 exists (serialize for reproducibility).
    if branch in ("paired", "separate"):
        paths["rna_done"] = str(INTERNAL / "artifacts" / "s1_rna_qc" / "rna_qc.h5ad")
    # For atac_only, S2 is the first modality stage after plan_review — demand
    # plan_review.md here so plan_review_propose is always pulled into the DAG.
    if branch == "atac_only":
        paths["plan_review_md"] = str(PRE_RUN / "summary" / "plan_review.md")
    return paths


rule s2_atac_qc_propose:
    input:
        unpack(_s2_propose_inputs)
    output:
        proposal = str(INTERNAL / "proposals" / "s2_atac_qc.yaml"),
        awaiting = str(INTERNAL / "proposals" / "s2_atac_qc.awaiting_approval"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        from executor import approval
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s2_atac_qc",
            "action": "TSS enrichment + n_fragments MAD via SnapATAC2 (no tile matrix here — S5 builds it)",
        }))
        approval.mark_awaiting(params.run_dir, "s2_atac_qc")


rule s2_atac_qc_execute:
    input:
        proposal         = str(INTERNAL / "proposals" / "s2_atac_qc.yaml"),
        approved         = str(INTERNAL / "checkpoints" / "s2_atac_qc.approved"),
        plan             = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done = str(INTERNAL / "checkpoints" / "plan_review.approved"),
    output:
        h5ad = str(INTERNAL / "artifacts" / "s2_atac_qc" / "atac_qc.h5ad"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s2_atac_qc"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s2_atac_qc", attempt),
        runtime=RUNTIME["s2_atac_qc"],
    run:
        import json
        from pathlib import Path
        from executor.stages import s2_atac_qc
        plan = json.loads(Path(input.plan).read_text())
        s2_atac_qc.run(params.run_dir, plan)
