"""Guard against the import-shadowing bug class.

A module-level ``import x`` re-imported *inside* a function turns ``x`` into a
function-local for that entire scope. Any use of ``x`` before the local import
line then raises ``UnboundLocalError`` at runtime — a deterministic crash that
unit tests on the happy path can miss if the local import sits on a rarely-hit
branch.

This exact bug took down S5 (`from .. import io as _io` re-imported mid-function
shadowed the module-level ``_io`` used earlier) and lurked latent in S3 (`sp`).
This test fails build-time if any ``executor/stages/*.py`` function re-imports a
name it already imports at module level, so the class cannot return.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

STAGES_DIR = Path(__file__).resolve().parent.parent / "executor" / "stages"


def _module_level_aliases(tree: ast.Module) -> set[str]:
    names: set[str] = set()
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                names.add(alias.asname or alias.name.split(".")[0])
    return names


def _shadowing_reimports(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(), str(path))
    modlevel = _module_level_aliases(tree)
    offences: list[str] = []
    for fn in ast.walk(tree):
        if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for node in ast.walk(fn):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                for alias in node.names:
                    name = alias.asname or alias.name.split(".")[0]
                    if name in modlevel:
                        offences.append(
                            f"{path.name}: function {fn.name!r} re-imports module-level "
                            f"name {name!r} at line {node.lineno} (shadows the module import)"
                        )
    return offences


@pytest.mark.parametrize("stage_file", sorted(STAGES_DIR.glob("*.py")), ids=lambda p: p.name)
def test_no_function_reimports_module_level_name(stage_file: Path) -> None:
    offences = _shadowing_reimports(stage_file)
    assert not offences, "Import-shadowing risk (UnboundLocalError class):\n" + "\n".join(offences)
