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


def test_state_model_records_supervisor_ownership_and_resume_source():
    state_model = (_root() / "contracts" / "state_model.md").read_text()
    assert "| `monitor.pid` | Processing `submit` / `supervisor-restart`" in state_model
    assert "| `latest_submission.json` | `execute-spec`" in state_model
    assert "resume-monitor reconstructs context" not in state_model
