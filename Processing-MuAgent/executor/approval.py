"""Checkpoint approval sentinels.

Each stage writes proposals/<stage>.awaiting_approval and checkpoints/<stage>.approved
as zero-byte or small YAML files. Snakemake depends on these sentinels.
"""
from __future__ import annotations

import getpass
import socket
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


def awaiting_path(run_dir: Path | str, stage: str) -> Path:
    from .run_paths import RunPaths
    return RunPaths(Path(run_dir)).awaiting_sentinel(stage)


def approved_path(run_dir: Path | str, stage: str) -> Path:
    from .run_paths import RunPaths
    return RunPaths(Path(run_dir)).approved_sentinel(stage)


def mark_awaiting(run_dir: Path | str, stage: str) -> Path:
    """Write the awaiting_approval sentinel and remove any prior approval.

    In auto-approve mode (env PMA_AUTO_APPROVE=1), we do NOT delete a pre-existing
    approval sentinel — the CLI has pre-seeded approvals so snakemake can execute
    the full DAG in one invocation without the propose rule invalidating them.
    """
    import os
    p = awaiting_path(run_dir, stage)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("")
    if os.environ.get("PMA_AUTO_APPROVE") == "1":
        return p
    approved = approved_path(run_dir, stage)
    if approved.exists():
        approved.unlink()
    return p


def approve(run_dir: Path | str, stage: str, *, actor: str | None = None,
            param_snapshot_hash: str | None = None, note: str = "") -> Path:
    p = approved_path(run_dir, stage)
    p.parent.mkdir(parents=True, exist_ok=True)
    record: dict[str, Any] = {
        "stage": stage,
        "actor": actor or f"{getpass.getuser()}@{socket.gethostname()}",
        "approved_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "param_snapshot_hash": param_snapshot_hash,
        "note": note,
    }
    with p.open("w") as f:
        yaml.safe_dump(record, f)
    awaiting = awaiting_path(run_dir, stage)
    if awaiting.exists():
        awaiting.unlink()
    return p


def is_approved(run_dir: Path | str, stage: str) -> bool:
    return approved_path(run_dir, stage).exists()
