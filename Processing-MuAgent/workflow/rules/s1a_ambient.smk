def _s1a_inputs(wildcards):
    # S1a only runs in branches that have RNA. For atac_only it produces an empty
    # pass-through. The dependency edge is s0's durable validation_report.json marker
    # (NOT the deletable raw RNA ingest h5ad, which S1a reads by path and which
    # _cleanup_qc_intermediates removes at the post_qc gate).
    return {
        "rna":   str(INTERNAL / "artifacts" / "s0_ingest" / "validation_report.json"),
        "plan":  str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        "plan_review_done": str(INTERNAL / "checkpoints" / "plan_review.approved"),
        "plan_review_md":   str(PLAN / f"plan_review_{RUN_DIR.name}.md"),
    }


rule s1a_ambient_propose:
    input:
        unpack(_s1a_inputs),
    output:
        proposal = str(INTERNAL / "proposals" / "s1a_ambient.yaml"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import json, yaml
        from pathlib import Path
        plan = json.loads(Path(input.plan).read_text())
        method_val = (plan.get("stages", {}).get("s1a_ambient", {})
                          .get("parameters", {}).get("method", {}).get("value", "auto"))
        # User revisions made at plan review take precedence over the frozen plan.
        params_path = Path(params.run_dir) / "internal" / "parameters.yaml"
        if params_path.exists():
            pdata = yaml.safe_load(params_path.read_text()) or {}
            override = pdata.get("s1a_ambient", {}).get("method", {}).get("value")
            if override:
                method_val = override
        _DESCRIPTIONS = {
            "auto":    "auto-select (SoupX if raw matrix present, else DecontX)",
            "decontx": "DecontX — filtered counts only (explicit)",
            "soupx":   "SoupX — raw + filtered counts (explicit)",
            "none":    "none — pass-through (ambient correction disabled)",
        }
        action = "Ambient RNA correction: " + _DESCRIPTIONS.get(
            str(method_val).lower(), str(method_val)
        )
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s1a_ambient",
            "action": action,
        }))


rule s1a_ambient_execute:
    input:
        plan     = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        # Durable S0 marker, not the deletable rna_ingest.h5ad (read by path in-stage).
        rna      = str(INTERNAL / "artifacts" / "s0_ingest" / "validation_report.json"),
    output:
        # summary.json is the SOLE declared output + durable stage-done marker (status
        # + the S1a->S1 edge key off it). rna_decontaminated.h5ad is an UNTRACKED
        # working file read by S1 by path and removed by _cleanup_qc_intermediates.
        summary = str(INTERNAL / "artifacts" / "s1a_ambient" / "summary.json"),
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
