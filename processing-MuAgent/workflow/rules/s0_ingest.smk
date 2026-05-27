rule s0_ingest_propose:
    input:
        approved_p1 = str(INTERNAL / "checkpoints" / "p1_context.approved"),
    output:
        proposal = str(INTERNAL / "proposals" / "s0_ingest.yaml"),
        awaiting = str(INTERNAL / "proposals" / "s0_ingest.awaiting_approval"),
    params:
        run_dir = str(RUN_DIR),
    run:
        import yaml
        from pathlib import Path
        from executor import approval, provenance
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
                        "auto-downgrade to 'separate'), persist optional translation parquet "
                        "for S2, handle metadata, hash inputs"),
        }))
        approval.mark_awaiting(params.run_dir, "s0_ingest")


rule s0_ingest_execute:
    input:
        proposal = str(INTERNAL / "proposals" / "s0_ingest.yaml"),
        approved = str(INTERNAL / "checkpoints" / "s0_ingest.approved"),
    output:
        report  = str(INTERNAL / "artifacts" / "s0_ingest" / "validation_report.json"),
        rna_h5  = str(INTERNAL / "artifacts" / "s0_ingest" / "rna_ingest.h5ad"),
    params:
        run_dir = str(RUN_DIR),
    run:
        from executor.stages import s0_ingest
        s0_ingest.run(params.run_dir, config)
