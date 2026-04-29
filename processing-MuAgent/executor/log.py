"""Structured append-only event log (log.jsonl)."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def log_event(run_dir: Path | str, event: dict[str, Any]) -> None:
    from .run_paths import RunPaths
    p = RunPaths(Path(run_dir)).log_jsonl
    p.parent.mkdir(parents=True, exist_ok=True)
    record = {"ts": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), **event}
    with p.open("a") as f:
        f.write(json.dumps(record) + "\n")
