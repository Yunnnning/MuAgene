"""Machine profile + manifest-driven env synthesis (the bootstrap-decoupling layer).

No real conda/registry needed — these exercise the on-disk machine.config contract and
the manifest-sourced environments section that lets provision/validate run with no
science site.config.
"""
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from execution_muagent import environment, machine
from execution_muagent.monitor import SiteConfig

_MANIFEST = (
    "schema_version: '1'\n"
    "cpu:\n  provider: lock\n  definition: workflow/envs/processing.yaml\n"
    "  lock: workflow/envs/processing.linux-64.lock\n"
    "  imports: workflow/envs/muagene.imports.txt\n"
    "gpu:\n  provider: container\n  definition: workflow/envs/muagene-gpu.def\n"
    "  imports: workflow/envs/muagene-gpu.imports.txt\n"
    "defaults:\n  gpu_image: ~/.muagene/images/muagene-gpu.sif\n"
)


class MachineConfigTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        envs = self.tmp / "workflow" / "envs"
        envs.mkdir(parents=True)
        (envs / "manifest.yaml").write_text(_MANIFEST)
        self.repo = self.tmp

    def test_write_load_roundtrip(self):
        p = self.tmp / "machine.config"
        cfg = machine.MachineConfig(processing_repo=str(self.repo), manager="mamba",
                                    singularity_module="singularityce/3.11.3",
                                    gpu_image_uri="docker://reg/muagene-gpu:25.04")
        machine.write_machine_config(cfg, p)
        loaded = machine.load_machine_config(p)
        self.assertEqual(loaded.manager, "mamba")
        self.assertEqual(loaded.gpu_image_uri, "docker://reg/muagene-gpu:25.04")
        self.assertEqual(loaded.singularity_module, "singularityce/3.11.3")

    def test_load_absent_returns_none(self):
        self.assertIsNone(machine.load_machine_config(self.tmp / "nope.config"))

    def test_default_environments_section_from_manifest(self):
        cfg = machine.MachineConfig(manager="mamba", singularity_module="m",
                                    gpu_image_uri="docker://reg/muagene-gpu:25.04")
        sec = machine.default_environments_section(self.repo, cfg)
        self.assertEqual(sec["cpu"]["provider"], "lock")
        self.assertEqual(sec["cpu"]["lock"], "workflow/envs/processing.linux-64.lock")
        self.assertEqual(sec["gpu"]["provider"], "container")
        self.assertEqual(sec["gpu"]["image_uri"], "docker://reg/muagene-gpu:25.04")
        self.assertEqual(sec["manager"], "mamba")

    def test_missing_manifest_raises(self):
        with self.assertRaises(FileNotFoundError):
            machine.load_env_manifest(self.tmp / "no_such_repo")

    def test_synthesize_site_config_resolves_specs(self):
        cfg = machine.MachineConfig(conda_env="muagene", gpu_conda_env="muagene-gpu",
                                    gpu_image_uri="docker://reg/muagene-gpu:25.04", manager="mamba")
        sc = machine.synthesize_site_config(self.repo, cfg)
        self.assertIsInstance(sc, SiteConfig)
        cpu = environment.resolve_env_spec(sc, self.repo, "cpu")
        self.assertEqual(cpu.provider, "lock")
        self.assertEqual(cpu.env_name, "muagene")
        self.assertTrue(str(cpu.lock).endswith("workflow/envs/processing.linux-64.lock"))
        gpu = environment.resolve_env_spec(sc, self.repo, "gpu")
        self.assertEqual(gpu.provider, "container")
        self.assertEqual(gpu.image_uri, "docker://reg/muagene-gpu:25.04")
        self.assertEqual(gpu.env_name, "muagene-gpu")

    def test_detect_machine_config_merges_overrides(self):
        with mock.patch.object(machine.capabilities, "probe_capabilities",
                               return_value={"manager": "conda", "container_runtime": "singularity",
                                             "scheduler": "slurm", "gpu_present": False}):
            cfg = machine.detect_machine_config(
                self.repo, manager="mamba", singularity_module="singularityce/3.11.3",
                gpu_image_uri="docker://reg/muagene-gpu:25.04")
        self.assertEqual(cfg.manager, "mamba")                  # explicit override wins
        self.assertEqual(cfg.container_runtime, "singularity")  # detected
        self.assertEqual(cfg.scheduler, "slurm")
        self.assertEqual(cfg.conda_env, "muagene")              # default identity
        self.assertEqual(cfg.gpu_conda_env, "muagene-gpu")


if __name__ == "__main__":
    unittest.main()
