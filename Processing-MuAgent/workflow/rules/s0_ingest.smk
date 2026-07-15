rule s0_ingest_propose:
    output:
        proposal = str(INTERNAL / "proposals" / "s0_ingest.yaml"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        from executor import provenance
        # Record config-supplied genome assembly as user-provenance
        provenance.set_param(
            str(INTERNAL / "parameters.yaml"),
            "ingest.genome_assembly",
            config.get("genome_assembly", "mm10"),
            source="user", confidence="high",
            rationale="Supplied by user in run.yaml",
        )
        _propose_keys = ("rna_path", "rna_raw_path", "atac_fragments_path",
                          "genome_assembly", "metadata_path",
                          "barcode_translation_path", "atac_peaks_path",
                          "cell_metadata_path")
        Path(output.proposal).write_text(yaml.safe_dump({
            "stage": "s0_ingest",
            "inputs": {k: config[k] for k in _propose_keys if k in config},
            "action": ("validate formats, detect pairing via the diagnostics ladder "
                        "(direct overlap -> suffix-normalized -> translation table -> "
                        "explicit branch confirmation), persist optional translation parquet "
                        "for S2, handle metadata, hash inputs"),
        }))


rule s0_ingest_execute:
    """Merged planning compute (single cluster job): load + validate + pair, assemble
    the preprocessing plan in-process, and run the QC threshold exploration on the
    in-memory matrices — emitting the data + figures plan_review consumes. Pulled in
    as a dependency when submit infers `plan_review_propose` (the planning-phase
    target); not a localrule because the ATAC fragment import + RNA QC are heavy.
    """
    input:
        context = str(INTERNAL / "artifacts" / "p1_context" / "context_extraction.json"),
    output:
        # validation_report.json is the durable S0 done-marker (also read post-gate by
        # S5). rna_ingest.h5ad and metadata_minimal.tsv are written as UNTRACKED working
        # files: read by S1a by path during the planning phase, then removed by
        # _cleanup_qc_intermediates at post_qc_review approval. Keeping them out of the
        # declared DAG means deleting them never triggers an S0 re-run — and rna_ingest.h5ad
        # is a pure cache (S1a reconstructs it via io.load_rna_ingest if a re-process needs it).
        report  = str(INTERNAL / "artifacts" / "s0_ingest" / "validation_report.json"),
        plan    = str(INTERNAL / "artifacts" / "p2_plan" / "preprocessing_plan.json"),
        explore = str(INTERNAL / "artifacts" / "qc_explore" / "qc_explore.json"),
    params:
        run_dir = str(RUN_DIR),
    threads: RESOURCES["s0_ingest"]["cpus"]
    resources:
        mem_mb=lambda wc, attempt: mem_mb_for("s0_ingest", attempt),
        runtime=RUNTIME["s0_ingest"],
    run:
        from executor.stages import s0_ingest
        from executor.cluster_exit import finalize_cluster_exit
        s0_ingest.run(params.run_dir, config)
        finalize_cluster_exit()
