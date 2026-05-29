"""Deliverables-directory bookkeeping.

Direct-write architecture: every logical artifact is written to its canonical
final location by the producing stage or helper. This module no longer
mirrors or symlinks files. It does two small things:

1. `ensure_scaffold(run_dir)` — create the deliverables/internal skeleton
   if missing. Idempotent.
2. `finalize(run_dir)` — sweep any legacy symlinks left by pre-refactor runs,
   then write a `layout.json` manifest listing the files actually present
   under `deliverables/`.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from .run_paths import RunPaths


def ensure_scaffold(run_dir: Path | str) -> None:
    RunPaths(Path(run_dir)).ensure()


def _clean_stale_symlinks(deliv_root: Path) -> int:
    """Remove any symlinks remaining from the pre-refactor mirrored layout.

    After the refactor every deliverable is a real file written directly by its
    producer. Any symlinks found here are leftovers from earlier pipeline
    revisions and are removed so fresh real-file writes land cleanly.
    """
    n_removed = 0
    if not deliv_root.exists():
        return 0
    for p in deliv_root.rglob("*"):
        if p.is_symlink():
            try:
                p.unlink()
                n_removed += 1
            except Exception:
                pass
    return n_removed


def _list_deliverables(deliv_root: Path) -> list[str]:
    """Return sorted relative paths of all real files under `deliverables/`."""
    if not deliv_root.exists():
        return []
    out: list[str] = []
    for p in deliv_root.rglob("*"):
        if p.is_file() and not p.is_symlink():
            out.append(str(p.relative_to(deliv_root)))
    out.sort()
    return out


def finalize(run_dir: Path | str) -> dict[str, Any]:
    """Sweep stale symlinks and write the `layout.json` manifest."""
    paths = RunPaths(Path(run_dir))
    paths.ensure()
    n_cleaned = _clean_stale_symlinks(paths.deliverables)
    deliverables_files = _list_deliverables(paths.deliverables)
    report: dict[str, Any] = {
        "layout_version": "2.0",
        "deliverables": deliverables_files,
        "stale_symlinks_cleaned": n_cleaned,
    }
    out = paths.layout_json
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    return report


# ---------------------------------------------------------------------------
# Backwards-compat shim — callers still invoke `layout.reorganise(run_dir)`;
# it now just calls `finalize`. No mirroring is performed.
# ---------------------------------------------------------------------------

def reorganise(run_dir: Path | str) -> dict[str, Any]:
    return finalize(run_dir)


def write_layout_report(run_dir: Path | str, report: dict[str, Any]) -> Path:
    """Kept for compatibility; layout.json is now written inside `finalize()`.
    This function just overwrites with the caller's report, preserving the API.
    """
    paths = RunPaths(Path(run_dir))
    out = paths.layout_json
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    return out
