rule manifest:
    input:
        s8 = str(INTERNAL / "artifacts" / "s8_umap" / "s8_done.txt"),
    output:
        manifest = str(RESULTS / "run_manifest.json"),
        qc_summary = str(RESULTS / "qc_summary.md"),
        notebook = str(RESULTS / "review_processed_h5mu.ipynb"),
        layout = str(RESULTS / "layout.json"),
    params:
        run_dir = str(RUN_DIR),
    run:
        from executor import manifest, provenance, layout, notebook_builder, qc_summary
        from pathlib import Path
        branch = provenance.current_branch(str(INTERNAL / "parameters.yaml"))
        cfg = dict(config)
        cfg["workflow_branch"] = branch
        cfg["run_id"] = cfg.get("run_id", Path(params.run_dir).name)
        # Write each deliverable directly to its canonical location. No mirroring.
        manifest.write_manifest(params.run_dir, cfg)
        qc_summary.write(params.run_dir)
        notebook_builder.write_review_notebook(params.run_dir)
        # Finalize: sweep any stale symlinks from pre-refactor runs + write
        # deliverables/results/layout.json manifest.
        layout.finalize(params.run_dir)
