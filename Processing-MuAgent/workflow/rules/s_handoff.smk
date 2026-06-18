rule s_handoff:
    """Terminal handoff bundle for Integration-MuAgent.

    Emits the per-sample post-QC deliverables — post_qc_<run>.h5mu (rna/atac mods,
    un-normalized post-doublet cells) + post_qc_manifest.json (the cross-package
    contract pointing at the RETAINED peaks BED + prepared ATAC fragments). Gated on
    post_qc_review approval; reads the S3 post-doublet artifacts (always written by
    S3, possibly empty on single-modality branches — s_handoff drops empty mods).

    Independently buildable terminal target (`run --target s_handoff`), ORTHOGONAL to
    S4–S8: it depends only on the S3 outputs + the post_qc_review gate, never on
    S4–S8. rule all requires both this bundle and run_manifest.json, so the existing
    end-to-end flow is untouched. When S4–S8 move to Integration-MuAgent, this becomes
    Preprocessing's terminus.
    """
    input:
        approved  = str(INTERNAL / "checkpoints" / "post_qc_review.approved"),
        plan      = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        rna_post  = str(INTERNAL / "artifacts" / "s3_doublets" / "rna_post_doublet.h5ad"),
        atac_post = str(INTERNAL / "artifacts" / "s3_doublets" / "atac_post_doublet.h5ad"),
    output:
        h5mu     = str(RESULTS / f"post_qc_{RUN_DIR.name}.h5mu"),
        manifest = str(RESULTS / "post_qc_manifest.json"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import json
        from pathlib import Path
        from executor.stages import s_handoff
        from executor import provenance
        plan = json.loads(Path(input.plan).read_text())
        branch = provenance.current_branch(str(INTERNAL / "parameters.yaml"))
        s_handoff.run(params.run_dir, plan, workflow_branch=branch)
