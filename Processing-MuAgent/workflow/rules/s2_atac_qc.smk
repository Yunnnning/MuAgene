def _s2_propose_inputs(wildcards):
    from executor import provenance
    branch = provenance.current_branch(str(INTERNAL / "parameters.yaml"))
    paths: dict = {
        "plan": str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
    }
    # Order-only dependency on S1 in branches where S1 exists (serialize for reproducibility).
    if branch in ("paired", "unpaired"):
        # Use qc_summary.json (not rna_qc.h5ad) — it survives post_qc_review cleanup.
        paths["rna_done"] = str(INTERNAL / "artifacts" / "s1_rna_qc" / "qc_summary.json")
    # For atac_only, S2 is the first modality stage after plan_review — demand
    # the run-scoped plan review md so plan_review_propose is always pulled into the DAG.
    if branch == "atac_only":
        paths["plan_review_md"] = str(PLAN / f"plan_review_{RUN_DIR.name}.md")
    return paths


rule s2_atac_qc_propose:
    input:
        unpack(_s2_propose_inputs)
    output:
        proposal = str(INTERNAL / "proposals" / "s2_atac_qc.yaml"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s2_atac_qc",
            "action": "TSS enrichment + n_fragments MAD via SnapATAC2 (no tile matrix here — S5 builds it)",
        }))


rule s2_atac_qc_execute:
    input:
        plan             = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done = str(INTERNAL / "checkpoints" / "plan_review.approved"),
    output:
        # qc_summary.json is the SOLE declared output and the durable stage-done
        # marker (status, reports, and the s3 dependency edge all key off it; the
        # h5ad read in qc_summary._atac_stat_rows is a guarded fallback).
        # atac_qc.h5ad is written by the stage as an UNTRACKED working file:
        # consumed only by s3_doublets (read by path) and removed by
        # _cleanup_qc_intermediates at post_qc_review approval. Keeping it out of
        # the declared DAG means deleting it never triggers a "Missing output
        # files" re-run of S2/S3.
        qc_summary = str(INTERNAL / "artifacts" / "s2_atac_qc" / "qc_summary.json"),
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
        from executor.cluster_exit import finalize_cluster_exit
        plan = json.loads(Path(input.plan).read_text())
        s2_atac_qc.run(params.run_dir, plan)
        finalize_cluster_exit()
