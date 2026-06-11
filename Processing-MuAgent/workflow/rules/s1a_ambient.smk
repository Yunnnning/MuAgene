def _s1a_inputs(wildcards):
    """S1a only runs in branches that have RNA. For atac_only it produces an
    empty pass-through; the dependency on s0's rna_ingest.h5ad is the only
    signal we need."""
    return {
        "rna":   str(INTERNAL / "artifacts" / "s0_ingest" / "rna_ingest.h5ad"),
        "plan":  str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        "plan_review_done": str(INTERNAL / "checkpoints" / "plan_review.approved"),
        "plan_review_md":   str(PLAN / "summary" / "plan_review.md"),
    }


rule s1a_ambient_propose:
    input:
        unpack(_s1a_inputs),
    output:
        proposal = str(INTERNAL / "proposals" / "s1a_ambient.yaml"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s1a_ambient",
            "action": (
                "Ambient RNA correction per plan (auto: SoupX if raw matrix present, "
                "else DecontX; none = pass-through; confirm method at plan review)"
            ),
        }))


rule s1a_ambient_execute:
    input:
        plan     = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        rna      = str(INTERNAL / "artifacts" / "s0_ingest" / "rna_ingest.h5ad"),
    output:
        h5ad = str(INTERNAL / "artifacts" / "s1a_ambient" / "rna_decontaminated.h5ad"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s1a_ambient"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s1a_ambient", attempt),
        runtime=RUNTIME["s1a_ambient"],
    run:
        import json
        from pathlib import Path
        from executor.stages import s1a_ambient
        from executor.cluster_exit import finalize_cluster_exit
        plan = json.loads(Path(input.plan).read_text())
        s1a_ambient.run(params.run_dir, plan)
        finalize_cluster_exit()
