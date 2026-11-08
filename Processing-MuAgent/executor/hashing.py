"""Input/parameter hashing and scientific tool-version collection."""
from __future__ import annotations

import hashlib
import json
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
