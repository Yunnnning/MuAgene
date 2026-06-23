"""Harness consistency tripwires (Execution side).

Guards the shared inter-agent finding-code registry: every finding code the
monitor/env commands emit must be registered in contracts/findings.yaml, and the
registry must cover the documented contract. Pure-Python, fast.
"""
from __future__ import annotations

import pathlib
import re

import yaml


def _root() -> pathlib.Path:
    # Execution-MuAgent/tests/<this> -> parents[2] == MuAgene repo root.
    return pathlib.Path(__file__).resolve().parents[2]


def _registry() -> set[str]:
    data = yaml.safe_load((_root() / "contracts" / "findings.yaml").read_text())
    return set(data["findings"])


# The finding-code contract as documented in Execution-MuAgent/README.md. The
# registry must cover all of these (a code removed from the registry, or a new
# documented code never registered, fails here).
DOCUMENTED_CODES = {
    "spec_validation_error", "submit_rejected_policy", "submit_rejected_transient",
    "scheduler_failed", "workflow_error_marker", "output_missing", "stall_confirmed",
    "filesystem_hang_suspected", "stall_suspected", "stall_recovered",
    "no_progress_files", "scheduler_completing", "scheduler_query_failed",
    "stage_output_verified", "workflow_complete",
    "env_missing", "env_stale", "env_stale_reprovision", "lock_stale_vs_yaml",
    "platform_unsupported", "gpu_image_unavailable", "gpu_import_needs_node",
    "import_failed", "provision_failed",
}


def test_registry_covers_documented_codes():
    missing = DOCUMENTED_CODES - _registry()
    assert not missing, f"codes documented but absent from contracts/findings.yaml: {sorted(missing)}"


def test_emitted_code_literals_are_registered():
    """Every `code: "x"` / `code="x"` literal in the package is registered."""
    reg = _registry()
    pkg = _root() / "Execution-MuAgent" / "execution_muagent"
    pat = re.compile(r"""code["']?\s*[:=]\s*["'](\w+)["']""")
    seen: set[str] = set()
    for py in pkg.glob("*.py"):
        seen |= set(pat.findall(py.read_text()))
    unregistered = seen - reg
    assert not unregistered, f"emitted finding codes not in registry: {sorted(unregistered)}"


def test_every_cli_command_is_documented():
    import importlib
    cli = importlib.import_module("execution_muagent.cli")
    grp = getattr(cli, "cli", None) or getattr(cli, "main", None)
    tools = (pathlib.Path(__file__).resolve().parents[1] / "agent" / "tools.md").read_text()
    missing = [c for c in grp.commands if c not in tools]
    assert not missing, f"Execution commands missing from agent/tools.md: {sorted(missing)}"
