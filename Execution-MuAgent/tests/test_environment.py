"""Unit tests for the env provisioning/validation + fingerprint contract.

Subprocess (build / create / import) is mocked, so these run anywhere with no conda,
container runtime, or GPU — the cross-machine-safe verification of the contract logic.
"""
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from execution_muagent import environment as env
from execution_muagent.monitor import SiteConfig


def _ok(stdout="", stderr=""):
    return subprocess.CompletedProcess([], 0, stdout, stderr)


def _fail(stderr="boom"):
    return subprocess.CompletedProcess([], 1, "", stderr)


class ResolveSpecTests(unittest.TestCase):
    def test_defaults_when_no_environments_section(self):
        # Back-compat: a site.config predating the contract still resolves.
        sc = SiteConfig(scheduler="slurm", conda_env="grn")
        cpu = env.resolve_env_spec(sc, "/repo", "cpu")
        gpu = env.resolve_env_spec(sc, "/repo", "gpu")
        self.assertEqual(cpu.provider, "yaml")   # no lock declared -> yaml fallback
        self.assertEqual(gpu.provider, "container")
        self.assertEqual(cpu.env_name, "grn")

    def test_reads_environments_section(self):
        sc = SiteConfig(scheduler="slurm", conda_env="muagene", gpu_conda_env="muagene-gpu",
                        environments={"singularity_module": "singularityce/3.11",
                                      "gpu": {"provider": "container",
                                              "definition": "workflow/envs/muagene-gpu.def",
                                              "image": "/imgs/x.sif",
                                              "image_uri": "docker://reg/muagene-gpu:25.04"}})
        gpu = env.resolve_env_spec(sc, "/repo", "gpu")
        self.assertEqual(gpu.provider, "container")
        self.assertEqual(str(gpu.image), "/imgs/x.sif")
        self.assertEqual(gpu.image_uri, "docker://reg/muagene-gpu:25.04")
        self.assertEqual(str(gpu.definition), "/repo/workflow/envs/muagene-gpu.def")
        self.assertEqual(gpu.singularity_module, "singularityce/3.11")


class FingerprintAndStatusTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self._state = self.tmp / "env_state.json"
        patcher = mock.patch.object(env, "_state_path", return_value=self._state)
        patcher.start()
        self.addCleanup(patcher.stop)

    def _container_spec(self, image_exists=True, image_uri="docker://reg/muagene-gpu:25.04"):
        sub = Path(tempfile.mkdtemp(dir=self.tmp))  # unique dir so specs don't clobber
        image = sub / "x.sif"
        if image_exists:
            image.write_text("fake-sif")
        return env.EnvSpec(device="gpu", provider="container", env_name="muagene-gpu",
                           definition=None, lock=None, image=image, imports=None,
                           image_uri=image_uri)

    def test_container_fingerprint_tracks_image_uri(self):
        # Pull-only: identity is the pinned image reference, not the .def content.
        s1 = self._container_spec(image_uri="docker://reg/muagene-gpu:25.04")
        s2 = self._container_spec(image_uri="docker://reg/muagene-gpu:25.06")
        self.assertTrue(env.compute_fingerprint(s1).startswith("uri:"))
        self.assertNotEqual(env.compute_fingerprint(s1), env.compute_fingerprint(s2))

    def test_lock_fingerprint_tracks_content(self):
        lock = self.tmp / "p.lock"
        lock.write_text("@EXPLICIT\nhttps://x/a.conda\n")
        spec = env.EnvSpec(device="cpu", provider="lock", env_name="muagene",
                           definition=None, lock=lock, image=None, imports=None)
        fp1 = env.compute_fingerprint(spec)
        self.assertTrue(fp1.startswith("sha256:"))
        lock.write_text("@EXPLICIT\nhttps://x/a.conda\nhttps://x/b.conda\n")
        self.assertNotEqual(fp1, env.compute_fingerprint(spec))

    def test_status_missing_then_ok_then_stale(self):
        spec = self._container_spec(image_exists=False)
        self.assertEqual(env.env_status(spec, manager=None), "missing")  # no image
        # provision (mock the PULL, create the image, record fingerprint)
        with mock.patch.object(env, "_bash_login",
                               side_effect=lambda *a, **k: (spec.image.write_text("sif"), _ok())[1]):
            res = env.provision_env(spec, SiteConfig(scheduler="slurm"), container_runtime="singularity")
        self.assertEqual(res["status"], "provisioned")
        self.assertEqual(res["action"], "pull_image")
        self.assertEqual(env.env_status(spec, manager=None), "ok")
        # republish a new image tag -> fingerprint drifts -> stale (image still present)
        spec.image_uri = "docker://reg/muagene-gpu:25.06"
        self.assertEqual(env.env_status(spec, manager=None), "stale")

    def test_provision_noop_when_ok(self):
        spec = self._container_spec(image_exists=True)
        env._record_provisioned(spec, env.compute_fingerprint(spec))
        with mock.patch.object(env, "_bash_login") as m:
            res = env.provision_env(spec, SiteConfig(scheduler="slurm"))
        self.assertEqual(res["action"], "noop")
        m.assert_not_called()

    def test_provision_fails_loud_without_image_uri(self):
        # Pull-only: a container with no image_uri cannot be provisioned (no local build).
        spec = self._container_spec(image_exists=False, image_uri=None)
        res = env.provision_env(spec, SiteConfig(scheduler="slurm"), container_runtime="singularity")
        self.assertEqual(res["status"], "failed")
        self.assertEqual(res.get("code"), "gpu_image_unavailable")


class ValidateAndReconcileTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        patcher = mock.patch.object(env, "_state_path", return_value=self.tmp / "state.json")
        patcher.start()
        self.addCleanup(patcher.stop)
        self.imports = self.tmp / "imp.txt"
        self.imports.write_text("rapids_singlecell\ncupy\n")
        self.image = self.tmp / "x.sif"
        self.image.write_text("sif")
        self.defn = self.tmp / "g.def"
        self.defn.write_text("Bootstrap: docker\n")
        self.spec = env.EnvSpec(device="gpu", provider="container", env_name="muagene-gpu",
                                definition=self.defn, lock=None, image=self.image, imports=self.imports,
                                image_uri="docker://reg/muagene-gpu:25.04")

    def test_validate_missing_env_is_error(self):
        self.image.unlink()
        res = env.validate_env(self.spec, SiteConfig(scheduler="slurm"))
        self.assertFalse(res["ok"])
        self.assertEqual(res["findings"][0]["code"], "env_missing")

    def test_validate_import_failure_is_error(self):
        with mock.patch.object(env, "_import_check", return_value=_fail("ModuleNotFoundError: cupy")):
            res = env.validate_env(self.spec, SiteConfig(scheduler="slurm"))
        self.assertFalse(res["ok"])
        self.assertTrue(any(f["code"] == "import_failed" for f in res["findings"]))

    def test_validate_cuda_unavailable_on_cpu_host_is_warning(self):
        with mock.patch.object(env, "_import_check",
                               return_value=_fail("cudaErrorNoDevice: no CUDA-capable device")), \
             mock.patch.object(env.capabilities, "gpu_present", return_value=False):
            res = env.validate_env(self.spec, SiteConfig(scheduler="slurm"))
        self.assertTrue(res["ok"])  # warning, not error
        self.assertTrue(any(f["code"] == "gpu_import_needs_node" for f in res["findings"]))

    def test_validate_ok_when_imports_succeed(self):
        env._record_provisioned(self.spec, env.compute_fingerprint(self.spec))
        with mock.patch.object(env, "_import_check", return_value=_ok()):
            res = env.validate_env(self.spec, SiteConfig(scheduler="slurm"))
        self.assertTrue(res["ok"])
        self.assertEqual(res["findings"], [])

    def test_reconcile_manual_policy_fails_loud_without_provisioning(self):
        self.image.unlink()  # missing
        sc = SiteConfig(scheduler="slurm", gpu_conda_env="muagene-gpu",
                        environments={"policy": "manual",
                                      "gpu": {"provider": "container",
                                              "definition": str(self.defn), "image": str(self.image),
                                              "imports": str(self.imports)}})
        with mock.patch.object(env, "_bash_login") as m:
            out = env.reconcile(sc, "/repo", "gpu")
        self.assertFalse(out["ok"])
        self.assertEqual(out["status"], "missing")
        m.assert_not_called()  # manual policy never builds
        self.assertTrue(any("provision-env" in f["message"] for f in out["findings"]))

    def test_reconcile_auto_policy_provisions_then_validates(self):
        self.image.unlink()  # missing -> auto should pull
        sc = SiteConfig(scheduler="slurm", gpu_conda_env="muagene-gpu",
                        environments={"policy": "auto",
                                      "gpu": {"provider": "container",
                                              "definition": str(self.defn), "image": str(self.image),
                                              "imports": str(self.imports),
                                              "image_uri": "docker://reg/muagene-gpu:25.04"}})

        def _pull(*a, **k):
            self.image.write_text("sif")
            return _ok()

        with mock.patch.object(env, "_bash_login", side_effect=_pull), \
             mock.patch.object(env, "_import_check", return_value=_ok()):
            out = env.reconcile(sc, "/repo", "gpu")
        self.assertTrue(out["ok"])
        self.assertEqual(out["provision"]["status"], "provisioned")
        self.assertEqual(out["provision"]["action"], "pull_image")


class LockPreflightTests(unittest.TestCase):
    """Linux-only + lock-freshness guards for the CPU (lock) provider."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        patcher = mock.patch.object(env, "_state_path", return_value=self.tmp / "state.json")
        patcher.start()
        self.addCleanup(patcher.stop)
        self.yaml = self.tmp / "processing.yaml"
        self.yaml.write_text("name: muagene\ndependencies: [scanpy]\n")
        self.lock = self.tmp / "processing.linux-64.lock"
        self._write_lock(self._yaml_hash(), "linux-64")  # fresh + matching

    def _yaml_hash(self):
        import hashlib
        return hashlib.sha256(self.yaml.read_bytes()).hexdigest()

    def _write_lock(self, source_hash, platform):
        self.lock.write_text(
            f"# source-sha256: {source_hash}\n# platform: {platform}\n@EXPLICIT\nhttps://x/a.conda\n")

    def _spec(self):
        return env.EnvSpec(device="cpu", provider="lock", env_name="muagene",
                           definition=self.yaml, lock=self.lock, image=None, imports=None)

    def _make_stale(self):
        # Edit the YAML so its content hash no longer matches the lock's recorded one
        # (git-stable: not mtime-based). Mirrors "edited yaml, forgot to regenerate".
        self.yaml.write_text("name: muagene\ndependencies: [scanpy, harmonypy]\n")

    def test_fresh_matching_lock_has_no_findings(self):
        with mock.patch.object(env, "_host_conda_subdir", return_value="linux-64"):
            self.assertEqual(env._lock_preflight_findings(self._spec()), [])

    def test_lock_without_marker_is_not_flagged(self):
        # A lock predating the source-sha256 convention can't be verified -> don't false-fire.
        self.lock.write_text("# platform: linux-64\n@EXPLICIT\nhttps://x/a.conda\n")
        self._make_stale()
        with mock.patch.object(env, "_host_conda_subdir", return_value="linux-64"):
            findings = env._lock_preflight_findings(self._spec())
        self.assertFalse(any(f["code"] == "lock_stale_vs_yaml" for f in findings))

    def test_lock_stale_vs_yaml_finding(self):
        self._make_stale()
        with mock.patch.object(env, "_host_conda_subdir", return_value="linux-64"):
            findings = env._lock_preflight_findings(self._spec())
        self.assertTrue(any(f["code"] == "lock_stale_vs_yaml" and f["severity"] == "error"
                            for f in findings))

    def test_platform_unsupported_on_non_linux(self):
        self._write_lock(self._yaml_hash(), "osx-arm64")
        with mock.patch.object(env, "_host_conda_subdir", return_value="linux-64"):
            res = env.validate_env(self._spec(), SiteConfig(scheduler="slurm"))
        self.assertFalse(res["ok"])
        self.assertEqual(res["findings"][0]["code"], "platform_unsupported")

    def test_validate_env_surfaces_lock_stale(self):
        self._make_stale()
        with mock.patch.object(env, "_host_conda_subdir", return_value="linux-64"), \
             mock.patch.object(env, "env_present", return_value=True), \
             mock.patch.object(env, "_import_check", return_value=_ok()):
            res = env.validate_env(self._spec(), SiteConfig(scheduler="slurm"))
        self.assertFalse(res["ok"])
        self.assertTrue(any(f["code"] == "lock_stale_vs_yaml" for f in res["findings"]))

    def test_reconcile_blocks_on_stale_lock_without_provisioning(self):
        self._make_stale()
        sc = SiteConfig(scheduler="slurm", conda_env="muagene",
                        environments={"policy": "auto", "cpu": {
                            "provider": "lock", "definition": str(self.yaml),
                            "lock": str(self.lock), "imports": None}})
        with mock.patch.object(env, "_host_conda_subdir", return_value="linux-64"), \
             mock.patch.object(env, "_run") as m_run, \
             mock.patch.object(env, "_bash_login") as m_bash:
            out = env.reconcile(sc, "/repo", "cpu")
        self.assertFalse(out["ok"])
        self.assertEqual(out["status"], "blocked")
        self.assertTrue(any(f["code"] == "lock_stale_vs_yaml" for f in out["findings"]))
        m_run.assert_not_called()
        m_bash.assert_not_called()


if __name__ == "__main__":
    unittest.main()
