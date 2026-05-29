"""sha256 of inputs, env fingerprint, param snapshot hashing."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path
from typing import Any


def sha256_file(path: Path | str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_json(obj: Any) -> str:
    return sha256_bytes(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode())


def env_fingerprint(tool_versions: dict[str, str]) -> str:
    payload = {"tool_versions": tool_versions, "python": os.sys.version.split()[0]}
    return sha256_json(payload)


def git_sha(repo_dir: Path | str) -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(repo_dir),
            capture_output=True,
            text=True,
            check=False,
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except FileNotFoundError:
        pass
    return "nogit"


def tool_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for mod in [
        "scanpy", "anndata", "muon", "mudata", "snapatac2", "scrublet",
        "leidenalg", "umap", "numpy", "scipy", "pandas",
    ]:
        try:
            m = __import__(mod)
            versions[mod] = getattr(m, "__version__", "unknown")
        except ImportError:
            versions[mod] = "MISSING"
    return versions
