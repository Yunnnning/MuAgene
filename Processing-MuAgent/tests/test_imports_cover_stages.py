"""Keep workflow/envs/*.imports.txt honest against what the stages actually import.

`validate-env` import-checks the env using these lists; if a stage starts importing a
library that is neither a declared dependency nor listed, `validate-env` would pass
while the stage fails at runtime. This static AST scan (no heavy imports) is the
tripwire. It runs anywhere — pure stdlib.
"""
import ast
import re
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ENVS = REPO / "workflow" / "envs"
SCAN = [*sorted((REPO / "executor" / "stages").glob("*.py")),
        REPO / "executor" / "io.py", REPO / "executor" / "compute.py"]

# pyproject distribution name -> import module name (only where they differ).
_DIST_TO_MODULE = {"pyyaml": "yaml", "umap-learn": "umap", "scikit-learn": "sklearn"}
# Heavy libs pulled in transitively (via scanpy/anndata), not direct pyproject deps.
_KNOWN_TRANSITIVE = {"h5py", "matplotlib"}
# Imported by child jobscripts (the Snakemake workflow), not by the stage modules.
_RUNTIME_ONLY = {"snakemake"}


def _imports_txt(name: str) -> set[str]:
    mods: set[str] = set()
    for line in (ENVS / name).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            mods.add(line)
    return mods


def _declared_modules() -> set[str]:
    text = (REPO / "pyproject.toml").read_text()
    block = re.search(r"dependencies\s*=\s*\[(.*?)\]", text, re.S)
    mods: set[str] = set()
    for raw in re.findall(r"['\"]([^'\"]+)['\"]", block.group(1) if block else ""):
        name = re.split(r"[<>=!~ \[]", raw, maxsplit=1)[0].strip().lower()
        mods.add(_DIST_TO_MODULE.get(name, name.replace("-", "_")))
    return mods


def _stage_imports() -> set[str]:
    mods: set[str] = set()
    for f in SCAN:
        if not f.exists():
            continue
        for node in ast.walk(ast.parse(f.read_text())):
            if isinstance(node, ast.Import):
                mods.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                mods.add(node.module.split(".")[0])  # absolute imports only
    return mods


class ImportsCoverStagesTests(unittest.TestCase):
    def test_no_dead_entries_in_imports_txt(self):
        used = _stage_imports() | _RUNTIME_ONLY
        for name in ("muagene.imports.txt", "muagene-gpu.imports.txt"):
            for mod in _imports_txt(name):
                self.assertIn(mod, used,
                    f"{name} lists {mod!r} but no stage imports it — stale entry or typo.")

    def test_gpu_distinguishers_present(self):
        # The whole point of the GPU list: a CPU-only env fails loud on these.
        self.assertTrue({"rapids_singlecell", "cupy"} <= _imports_txt("muagene-gpu.imports.txt"))

    def test_every_stage_import_is_declared_or_listed(self):
        allowed = (_declared_modules() | _imports_txt("muagene.imports.txt")
                   | _imports_txt("muagene-gpu.imports.txt") | _KNOWN_TRANSITIVE
                   | set(sys.stdlib_module_names))
        offenders = {m for m in _stage_imports() if m not in allowed and not m.startswith("_")}
        self.assertEqual(offenders, set(),
            f"stage modules not in pyproject deps / imports.txt / known-transitive: {sorted(offenders)}. "
            "Declare the dependency (and add to imports.txt if it is env-distinguishing).")


if __name__ == "__main__":
    unittest.main()
