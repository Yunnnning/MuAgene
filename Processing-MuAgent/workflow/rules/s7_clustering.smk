rule s7_clustering_propose:
    """Run the resolution sweep on both modalities; write recommendation + sweep table + figures.
    Final cluster labels are NOT assigned here — that happens in s7_clustering_execute after
    user approves (or revises) s7_clustering.rna.resolution / s7_clustering.atac.resolution
    in parameters.yaml.
    """
    input:
        rna_neighbors = str(INTERNAL / "artifacts" / "s6_neighbors" / "rna_neighbors.h5ad"),
        plan       = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
    output:
        proposal = str(INTERNAL / "proposals" / "s7_clustering.yaml"),
        awaiting = str(INTERNAL / "proposals" / "s7_clustering.awaiting_approval"),
        sweep    = str(INTERNAL / "artifacts" / "s7_clustering" / "sweep.parquet"),
        # resolution_summary.md is user-facing → deliverables/summary/
        summary  = str(CHECKPOINTS / "resolution_review" / "resolution_summary.md"),
        # Review notebook + static HTML deliverables (built statically from the
        # sweep artifacts; safe to declare as outputs since they are always written).
        notebook = str(CHECKPOINTS / "resolution_review" / "resolution_review.ipynb"),
        html     = str(CHECKPOINTS / "resolution_review" / "resolution_review.html"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import json
        import yaml
        from pathlib import Path
        from executor import approval, provenance
        from executor.stages import s7_clustering, s7_notebook
        plan = json.loads(Path(input.plan).read_text())
        result = s7_clustering.propose(params.run_dir, plan)
        # Build the resolution-review notebook + static HTML; this is the primary
        # deliverable users open during the S7 pause in headless / HPC runs.
        s7_notebook.build_and_render(params.run_dir)
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s7_clustering",
            "action": "resolution sweep complete; review resolution_review.html (or .ipynb), "
                      "then approve or revise s7_clustering.rna.resolution / "
                      "s7_clustering.atac.resolution",
            "recommended": result,
            "review_artifacts": {
                "html": str(Path(output.html).resolve()),
                "notebook": str(Path(output.notebook).resolve()),
                "markdown": str(Path(output.summary).resolve()),
            },
        }))
        approval.mark_awaiting(params.run_dir, "s7_clustering")


rule s7_clustering_execute:
    """Read approved resolutions from parameters.yaml and assign final cluster labels.
    Blocks until checkpoints/s7_clustering.approved exists.
    """
    input:
        proposal         = str(INTERNAL / "proposals" / "s7_clustering.yaml"),
        sweep            = str(INTERNAL / "artifacts" / "s7_clustering" / "sweep.parquet"),
        approved         = str(INTERNAL / "checkpoints" / "s7_clustering.approved"),
        plan             = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done = str(INTERNAL / "checkpoints" / "plan_review.approved"),
        rna_neighbors    = str(INTERNAL / "artifacts" / "s6_neighbors" / "rna_neighbors.h5ad"),
    output:
        rna_clustered = str(INTERNAL / "artifacts" / "s7_clustering" / "rna_clustered.h5ad"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s7_clustering"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s7_clustering", attempt),
        runtime=RUNTIME["s7_clustering"],
    run:
        import json
        from pathlib import Path
        from executor.stages import s7_clustering
        from executor.cluster_exit import finalize_cluster_exit
        plan = json.loads(Path(input.plan).read_text())
        s7_clustering.execute(params.run_dir, plan)
        finalize_cluster_exit()
