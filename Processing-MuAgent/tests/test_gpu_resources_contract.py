"""Contract test for the GPU-resource growth path (workflow/resources.smk).

Snakemake `default-resources` cannot key off the rule name, so a stage that is in
`_GPU_CAPABLE` but whose `<stage>_execute` rule omits
`gpu=RESOURCES["<stage>"]["gpu"]` would silently fall back to the profile default
(gpu=0) and request NO GPU. This test enforces the 3-edit growth contract documented
in resources.smk: every _GPU_CAPABLE stage must declare its own `gpu=` resource.

The set is empty today (preprocessing is CPU-only), so the contract check is
vacuously green; the second test verifies the checker itself so it fails loud the
moment the integration subagent adds a GPU stage without the declaration.
"""
import importlib.machinery
import importlib.util
import re
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RESOURCES_SMK = REPO / "workflow" / "resources.smk"
RULES_DIR = REPO / "workflow" / "rules"


def _load_gpu_capable():
    """Import _GPU_CAPABLE from the .smk file (it is plain, Snakemake-free Python)."""
    loader = importlib.machinery.SourceFileLoader("pma_resources_smk", str(RESOURCES_SMK))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod._GPU_CAPABLE


def _execute_rule_block(stage, text):
    """Text of `rule <stage>_execute:` up to the next top-level `rule `/EOF (or None)."""
    m = re.search(rf"^rule {re.escape(stage)}_execute:.*?(?=^rule \w|\Z)",
                  text, re.MULTILINE | re.DOTALL)
    return m.group(0) if m else None


def _declares_gpu(block):
    """True if the rule block declares a `gpu=RESOURCES[...]` resource."""
    return bool(re.search(r"gpu\s*=\s*RESOURCES\[", block))


class GpuResourceContractTests(unittest.TestCase):
    def test_every_gpu_capable_stage_declares_gpu_resource(self):
        gpu_capable = _load_gpu_capable()
        all_rules = "\n".join(p.read_text() for p in sorted(RULES_DIR.glob("*.smk")))
        for stage in sorted(gpu_capable):
            block = _execute_rule_block(stage, all_rules)
            self.assertIsNotNone(
                block, f"no `rule {stage}_execute:` found for _GPU_CAPABLE stage {stage!r}")
            self.assertTrue(
                _declares_gpu(block),
                f"_GPU_CAPABLE stage {stage!r} must declare "
                f'gpu=RESOURCES["{stage}"]["gpu"] in its {stage}_execute resources: block '
                "(else Snakemake's default gpu=0 applies and it silently requests 0 GPUs).")

    def test_checker_detects_present_and_missing_declaration(self):
        # Guards the safety net itself while _GPU_CAPABLE is empty.
        present = (
            'rule foo_execute:\n'
            '    resources:\n'
            '        mem_mb=1000,\n'
            '        gpu=RESOURCES["foo"]["gpu"],\n'
        )
        missing = (
            'rule foo_execute:\n'
            '    resources:\n'
            '        mem_mb=1000,\n'
        )
        self.assertTrue(_declares_gpu(_execute_rule_block("foo", present)))
        self.assertFalse(_declares_gpu(_execute_rule_block("foo", missing)))


if __name__ == "__main__":
    unittest.main()
