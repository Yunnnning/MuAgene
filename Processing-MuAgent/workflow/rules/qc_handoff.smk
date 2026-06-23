rule qc_handoff:
    """Terminal post-QC handoff bundle for Integration-MuAgent.

    Emits the per-sample post-QC deliverables — post_qc_<run>.h5mu (rna/atac mods,
    un-normalized post-doublet cells) + post_qc_manifest.json (the cross-package
    contract pointing at the RETAINED peaks BED + prepared ATAC fragments). Gated on
    post_qc_review approval; reads the S3 post-doublet artifacts (always written by
    S3) by path. On single-modality branches the absent side is a degenerate empty
    placeholder and is dropped; a modality the branch *expects* but cannot load is a
    hard error (no silent partial bundle). As its last step it DELETES the transient
    post-doublet h5ads — the h5mu is now the canonical post-QC store and S4/S5 read it;
    calls.parquet (the declared dep), joint_barcodes, overlap_summary, peaks and
    fragments are kept.

    NOT a localrule: the ATAC side is a SnapATAC2 (Blosc-compressed) matrix that must
    be read via snap.read and re-encoded portably, which is too heavy for the
    login/head node. On HPC this runs as a SLURM job; under `run` (local) it executes
    inline.

    Independently buildable (`run --target qc_handoff`) and depends only on the durable
    S3 marker (calls.parquet) + the post_qc_review gate — never on S4–S8. It is now
    UPSTREAM of S4/S5, which read the post-QC h5mu it writes (the internal post-doublet
    h5ads it deletes are gone). Running it at QC approval produces the h5mu early; the
    finish-batch `all` then finds it up-to-date and skips it. rule all requires both
    this bundle and run_manifest.json. When S4–S8 move to Integration-MuAgent, this
    becomes Preprocessing's terminus.
    """
    input:
        approved  = str(INTERNAL / "checkpoints" / "post_qc_review.approved"),
        plan      = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        # Durable S3 marker — NOT the post-doublet h5ads, which qc_handoff reads by
        # path and then deletes. calls.parquet survives that deletion and carries the
        # S3 -> qc_handoff edge (re-running S3 regenerates it, re-triggering qc_handoff).
        calls     = str(INTERNAL / "artifacts" / "s3_doublets" / "calls.parquet"),
    output:
        h5mu     = str(QC / f"post_qc_{RUN_DIR.name}.h5mu"),
        manifest = str(QC / "post_qc_manifest.json"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["qc_handoff"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("qc_handoff", attempt),
        runtime=RUNTIME["qc_handoff"],
    run:
        import json
        from pathlib import Path
        from executor.stages import qc_handoff
        from executor import provenance
        from executor.cluster_exit import finalize_cluster_exit
        plan = json.loads(Path(input.plan).read_text())
        branch = provenance.current_branch(str(INTERNAL / "parameters.yaml"))
        qc_handoff.run(params.run_dir, plan, workflow_branch=branch)
        finalize_cluster_exit()
