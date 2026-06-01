rule s8_umap_propose:
    input:
        rna_clustered = str(INTERNAL / "artifacts" / "s7_clustering" / "rna_clustered.h5ad"),
    output:
        proposal = str(INTERNAL / "proposals" / "s8_umap.yaml"),
        awaiting = str(INTERNAL / "proposals" / "s8_umap.awaiting_approval"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        from executor import approval
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s8_umap",
            "action": "per-modality UMAP; final h5mu (paired) or two h5ad (separate). HARD STOP.",
        }))
        approval.mark_awaiting(params.run_dir, "s8_umap")


def _s8_outputs(wildcards):
    from executor import provenance
    branch = provenance.get_value(str(INTERNAL / "parameters.yaml"),
                                  "plan.workflow_branch", "paired")
    if branch == "paired":
        return [str(POST_RUN / "processed.h5mu")]
    return [
        str(POST_RUN / "rna_processed.h5ad"),
        str(POST_RUN / "atac_processed.h5ad"),
    ]


rule s8_umap_execute:
    input:
        proposal         = str(INTERNAL / "proposals" / "s8_umap.yaml"),
        approved         = str(INTERNAL / "checkpoints" / "s8_umap.approved"),
        plan             = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done = str(INTERNAL / "checkpoints" / "plan_review.approved"),
        rna_clustered    = str(INTERNAL / "artifacts" / "s7_clustering" / "rna_clustered.h5ad"),
    output:
        # We always produce a sentinel file; branch-specific outputs are also written.
        sentinel = str(INTERNAL / "artifacts" / "s8_umap" / "s8_done.txt"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s8_umap"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s8_umap", attempt),
        runtime=RUNTIME["s8_umap"],
    run:
        import json
        from pathlib import Path
        from executor.stages import s8_umap
        from executor import provenance
        plan = json.loads(Path(input.plan).read_text())
        branch = provenance.get_value(str(INTERNAL / "parameters.yaml"),
                                      "plan.workflow_branch", "paired")
        result = s8_umap.run(params.run_dir, plan, workflow_branch=branch)
        Path(output.sentinel).write_text(str(result))
