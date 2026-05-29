#!/usr/bin/env python3
"""Apply the NFS-hang fix to snakemake-executor-plugin-slurm.

Run once after creating / updating the conda environment:

    python scripts/apply_nfs_patch.py

Background
----------
snakemake-executor-plugin-slurm v2.x submits each child SLURM job with
``--executor slurm-jobstep``.  The jobstep plugin sets ``ExecMode.REMOTE``
on the child snakemake process, which triggers ``dag.store_storage_outputs()``
after every rule finishes.  On NFS-mounted workdirs that function enters an
async event loop that blocks indefinitely — even when no files need
transferring — because every ``await`` yields to a background coroutine that
is stuck on a slow NFS stat.  The result: child SLURM jobs whose actual rule
completes in seconds are killed by the walltime limit.

Fix: change the child-job executor argument from ``--executor slurm-jobstep``
to ``--executor local``.  With ``ExecMode.DEFAULT`` (local executor) the child
snakemake process skips ``store_storage_outputs()`` entirely.  The rule still
runs on the allocated compute node; only the srun-based jobstep wrapper is
removed.  For NFS-backed clusters where all data is already on a shared
filesystem this is the correct execution model.
"""
from __future__ import annotations

import sys
import importlib.util
from pathlib import Path


TARGET_OLD = 'general_args = "--executor slurm-jobstep --jobs 1"'
TARGET_NEW = 'return "--executor local --cores 4"  # NFS patch: skip slurm-jobstep'
MARKER = "# NFS patch"


def find_plugin_file() -> Path:
    spec = importlib.util.find_spec("snakemake_executor_plugin_slurm")
    if spec is None or spec.origin is None:
        sys.exit(
            "ERROR: snakemake_executor_plugin_slurm not found. "
            "Activate the grn env and re-run."
        )
    return Path(spec.origin)


def patch(path: Path) -> None:
    text = path.read_text()

    if MARKER in text:
        print(f"Already patched: {path}")
        return

    if TARGET_OLD not in text:
        sys.exit(
            f"ERROR: expected patch target not found in {path}.\n"
            f"The plugin may have changed — check that you have version 2.7.0 "
            f"(pip show snakemake-executor-plugin-slurm) and review the diff "
            f"against the TARGET_OLD string in this script."
        )

    patched = text.replace(TARGET_OLD, TARGET_NEW, 1)
    path.write_text(patched)
    print(f"Patched: {path}")
    print("  slurm-jobstep → local executor for NFS-safe child job execution.")


if __name__ == "__main__":
    plugin_file = find_plugin_file()
    patch(plugin_file)
