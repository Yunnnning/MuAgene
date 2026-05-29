def _s6_inputs(wildcards):
    from executor import provenance
    branch = provenance.current_branch(str(INTERNAL / "parameters.yaml"))
    paths: dict = {}
    if branch in ("paired", "separate", "rna_only"):
        paths["rna_norm"] = str(INTERNAL / "artifacts" / "s4_rna_norm" / "rna_norm.h5ad")
    if branch in ("paired", "separate", "atac_only"):
        paths["atac_lsi"] = str(INTERNAL / "artifacts" / "s5_atac_lsi" / "lsi_summary.json")
    return paths


rule s6_dimred_propose:
    input:
        unpack(_s6_inputs)
    output:
        proposal = str(INTERNAL / "proposals" / "s6_dimred.yaml"),
        awaiting = str(INTERNAL / "proposals" / "s6_dimred.awaiting_approval"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        from executor import approval
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s6_dimred",
            "action": "RNA PCA + neighbors; ATAC neighbors on spectral/LSI components 2..N",
        }))
        approval.mark_awaiting(params.run_dir, "s6_dimred")


rule s6_dimred_execute:
    input:
        proposal         = str(INTERNAL / "proposals" / "s6_dimred.yaml"),
        approved         = str(INTERNAL / "checkpoints" / "s6_dimred.approved"),
        plan             = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        plan_review_done = str(INTERNAL / "checkpoints" / "plan_review.approved"),
    output:
        rna_h5ad = str(INTERNAL / "artifacts" / "s6_dimred" / "rna_dimred.h5ad"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s6_dimred"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s6_dimred", attempt),
        runtime=RUNTIME["s6_dimred"],
    run:
        import json
        from pathlib import Path
        from executor.stages import s6_dimred
        plan = json.loads(Path(input.plan).read_text())
        s6_dimred.run(params.run_dir, plan)
