#!/usr/bin/env python3
"""Apply the NFS-hang fix to snakemake-executor-plugin-slurm.

Run once after creating / updating the conda environment:

    python scripts/apply_nfs_patch.py

Background
----------
snakemake-executor-plugin-slurm v2.x submits each child SLURM job with
``--executor slurm-jobstep``.  The jobstep executor sets ``ExecMode.REMOTE``
on the child snakemake process via two mechanisms:

  1. ``additional_general_args()`` returns ``--executor slurm-jobstep --jobs 1``,
     which resolves to the jobstep plugin.
  2. ``get_exec_mode()`` (inherited from RemoteExecutor) returns
     ``ExecMode.REMOTE``, injecting ``--mode remote`` into every child job
     command line.

With ``ExecMode.REMOTE``, ``workflow.remote_exec`` is True, which causes
Snakemake to call ``dag.store_storage_outputs()`` after every rule finishes.
On NFS-mounted workdirs that function enters an async event loop that blocks
indefinitely — even when no files need transferring — because every ``await``
yields to a background coroutine stuck on a slow NFS stat.  Child SLURM jobs
whose actual rule completes in seconds are killed by the walltime limit.

Fix (two-part)
--------------
1. ``additional_general_args()``: return ``--executor local --cores 4``
   instead of ``--executor slurm-jobstep --jobs 1``.
2. ``get_exec_mode()``: override the inherited method to return
   ``ExecMode.DEFAULT`` instead of ``ExecMode.REMOTE``.  Without this,
   ``--mode remote`` is still injected into child jobs regardless of the
   executor change, keeping ``remote_exec=True`` and still triggering
   ``store_storage_outputs()``.

Together these ensure child jobs run with ``ExecMode.DEFAULT``, so
``store_storage_outputs()`` is never called.  Rules run directly on the
allocated compute node; only the srun-based jobstep wrapper is removed.
This is correct for NFS-backed clusters where all data is already on a
shared filesystem.
"""
from __future__ import annotations

import sys
import importlib.util
from pathlib import Path


# ---------------------------------------------------------------------------
# Patch 1 – replace additional_general_args() body
#   Handles three cases:
#     a) original unpatched code (TARGET_P1_ORIG)
#     b) manually patched without marker (TARGET_P1_BARE)
#     c) already fully patched (TARGET_P1_MARKED) — skip
# ---------------------------------------------------------------------------
TARGET_P1_ORIG   = 'general_args = "--executor slurm-jobstep --jobs 1"'
TARGET_P1_BARE   = '        return "--executor local --cores 4"'
TARGET_P1_MARKED = '        return "--executor local --cores 4"  # NFS patch: skip slurm-jobstep'

# ---------------------------------------------------------------------------
# Patch 2 – add ExecMode import
# ---------------------------------------------------------------------------
TARGET_P2_OLD = (
    "from snakemake_interface_executor_plugins.settings import (\n"
    "    ExecutorSettingsBase,\n"
    "    CommonSettings,\n"
    ")"
)
TARGET_P2_NEW = (
    "from snakemake_interface_executor_plugins.settings import (\n"
    "    ExecutorSettingsBase,\n"
    "    CommonSettings,\n"
    "    ExecMode,  # NFS patch: needed for get_exec_mode override\n"
    ")"
)
MARKER_P2 = "NFS patch: needed for get_exec_mode override"

# ---------------------------------------------------------------------------
# Patch 3 – add get_exec_mode() override after additional_general_args()
# ---------------------------------------------------------------------------
# The correct ExecMode for child jobs is SUBPROCESS:
#   - Skips store_storage_outputs() (remote_exec=False)
#   - Skips cleanup_workdir() (only called for DEFAULT)
#   - Disables file logging in child (enable_file_logging=False) — avoids NFS log writes
#   - Sets store_in_storage=False in per-job postprocess (same as REMOTE)
#   - Skips persistence metadata recording (only for DEFAULT)
# DEFAULT was tried first but still hung because cleanup_workdir() does NFS scandir/stat
# on every job's input/output files after rule completion.
TARGET_P3_OLD = (
    '        return "--executor local --cores 4"  # NFS patch: skip slurm-jobstep\n'
    "\n"
    "    def run_jobs("
)
# Handle migration from a previous wrong patch (DEFAULT → SUBPROCESS):
TARGET_P3_OLD_WRONG = (
    '        return "--executor local --cores 4"  # NFS patch: skip slurm-jobstep\n'
    "\n"
    "    def get_exec_mode(self):"
    "  # NFS patch: DEFAULT so child jobs don't set remote_exec=True\n"
    "        return ExecMode.DEFAULT\n"
    "\n"
    "    def run_jobs("
)
TARGET_P3_NEW = (
    '        return "--executor local --cores 4"  # NFS patch: skip slurm-jobstep\n'
    "\n"
    "    def get_exec_mode(self):"
    "  # NFS patch: SUBPROCESS skips store_storage_outputs, cleanup_workdir, NFS log file\n"
    "        return ExecMode.SUBPROCESS\n"
    "\n"
    "    def run_jobs("
)
MARKER_P3 = "NFS patch: SUBPROCESS skips store_storage_outputs, cleanup_workdir, NFS log file"
MARKER_P3_WRONG = "NFS patch: DEFAULT so child jobs don't set remote_exec=True"


def find_plugin_file() -> Path:
    spec = importlib.util.find_spec("snakemake_executor_plugin_slurm")
    if spec is None or spec.origin is None:
        sys.exit(
            "ERROR: snakemake_executor_plugin_slurm not found. "
            "Activate the grn env and re-run."
        )
    return Path(spec.origin)


def apply_patch(path: Path) -> None:
    text = path.read_text()
    changed = False

    # --- Patch 1 ---
    if TARGET_P1_MARKED in text:
        print("  P1 (additional_general_args): already patched")
    elif TARGET_P1_BARE in text:
        # Manually patched without marker — add marker
        text = text.replace(TARGET_P1_BARE, TARGET_P1_MARKED, 1)
        print("  P1 (additional_general_args): marker added to existing fix")
        changed = True
    elif TARGET_P1_ORIG in text:
        # Original unpatched code
        text = text.replace(TARGET_P1_ORIG, TARGET_P1_MARKED, 1)
        print("  P1 (additional_general_args): applied (slurm-jobstep → local)")
        changed = True
    else:
        sys.exit(
            f"ERROR: P1 patch target not found in {path}.\n"
            "The plugin may have changed — check snakemake-executor-plugin-slurm "
            "version (expected 2.7.0) and review TARGET_P1_* strings in this script."
        )

    # --- Patch 2 (ExecMode import) ---
    if MARKER_P2 in text:
        print("  P2 (ExecMode import): already patched")
    elif TARGET_P2_OLD in text:
        text = text.replace(TARGET_P2_OLD, TARGET_P2_NEW, 1)
        print("  P2 (ExecMode import): added ExecMode to settings import")
        changed = True
    else:
        sys.exit(
            f"ERROR: P2 patch target (settings import block) not found in {path}.\n"
            "Review TARGET_P2_OLD string in this script."
        )

    # --- Patch 3 (get_exec_mode override) ---
    if MARKER_P3 in text:
        print("  P3 (get_exec_mode): already patched (SUBPROCESS)")
    elif MARKER_P3_WRONG in text:
        # Previous patch used DEFAULT — migrate to SUBPROCESS
        text = text.replace(TARGET_P3_OLD_WRONG, TARGET_P3_NEW, 1)
        print("  P3 (get_exec_mode): migrated DEFAULT → SUBPROCESS")
        changed = True
    elif TARGET_P3_OLD in text:
        text = text.replace(TARGET_P3_OLD, TARGET_P3_NEW, 1)
        print("  P3 (get_exec_mode): added override returning ExecMode.SUBPROCESS")
        changed = True
    else:
        sys.exit(
            f"ERROR: P3 patch target (after additional_general_args return line) "
            f"not found in {path}.\n"
            "Review TARGET_P3_OLD / TARGET_P3_OLD_WRONG strings in this script."
        )

    if changed:
        path.write_text(text)
        print(f"\nPatched: {path}")
    else:
        print(f"\nNothing to do — all patches already applied: {path}")


if __name__ == "__main__":
    plugin_file = find_plugin_file()
    print(f"Plugin: {plugin_file}")
    apply_patch(plugin_file)
