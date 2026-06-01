rule p1_context_propose:
    output:
        proposal = str(INTERNAL / "proposals" / "p1_context.yaml"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        from executor import context as ctx
        # Canonical user-facing location for the Biological Context Report
        default_report_path = str(PRE_RUN / "config" / "biological_context.md")
        report_path = config.get("biological_context_path") or default_report_path
        ctx.write_template(report_path)
        Path(output.proposal).parent.mkdir(parents=True, exist_ok=True)
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "p1_context",
            "report_path": report_path,
            "action": "run executor.context.extract_context after user fills the report",
        }))


rule p1_context_execute:
    output:
        # Machine-readable extraction stays internal; human-readable summary is
        # a pre-run deliverable (reviewed before plan approval).
        extraction = str(INTERNAL / "artifacts" / "p1_context" / "context_extraction.json"),
        summary    = str(PRE_RUN / "summary" / "context_summary.md"),
    params:
        run_dir = str(RUN_DIR),
    run:
        from pathlib import Path
        from executor import context as ctx
        default_report_path = str(PRE_RUN / "config" / "biological_context.md")
        report_path = config.get("biological_context_path") or default_report_path
        payload = ctx.extract_context(report_path, params.run_dir,
                                      file_input_signals={})
        ctx.write_context_extraction(params.run_dir, payload)
        Path(output.summary).parent.mkdir(parents=True, exist_ok=True)
        Path(output.summary).write_text(ctx.render_summary(payload))
