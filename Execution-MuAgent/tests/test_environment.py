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


class CreateFromLockTests(unittest.TestCase):
    """Tests for the create-vs-update routing and broken-prefix cleanup in _create_from_lock."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.lock = self.tmp / "processing.linux-64.lock"
        self.lock.write_text("@EXPLICIT\nhttps://x/a.conda\n")
        self.spec = env.EnvSpec(device="cpu", provider="lock", env_name="grn",
                                definition=None, lock=self.lock, image=None, imports=None)

    def test_creates_when_env_absent(self):
        with mock.patch.object(env, "_named_env_usable", return_value=False), \
             mock.patch.object(env, "_clean_broken_env_prefix") as m_clean, \
             mock.patch.object(env, "_run", return_value=_ok()) as m_run:
            result = env._create_from_lock(self.spec, "micromamba")
        self.assertEqual(result["action"], "create_from_lock")
        self.assertIn("create", m_run.call_args[0][0])
        m_clean.assert_called_once_with("micromamba", "grn")

    def test_updates_in_place_when_env_usable(self):
        with mock.patch.object(env, "_named_env_usable", return_value=True), \
             mock.patch.object(env, "_run", return_value=_ok()) as m_run:
            result = env._create_from_lock(self.spec, "micromamba")
        self.assertEqual(result["action"], "update_from_lock")
        cmd = m_run.call_args[0][0]
        self.assertIn("install", cmd)
        self.assertNotIn("create", cmd)

    def test_recreates_broken_listed_env(self):
        # Listed by `env list` but the prefix is broken (no conda-meta/history): must route to
        # create (after cleanup), NOT update — updating a half-built env can't repair it.
        envs_dir = self.tmp / "envs"
        broken = envs_dir / "grn"
        (broken / "conda-meta").mkdir(parents=True)  # no history -> broken
        with mock.patch.object(env, "_conda_env_present", return_value=True), \
             mock.patch.object(env, "_conda_envs_dirs", return_value=[envs_dir]), \
             mock.patch.object(env, "_run", return_value=_ok()) as m_run:
            result = env._create_from_lock(self.spec, "micromamba")
        self.assertEqual(result["action"], "create_from_lock")
        self.assertIn("create", m_run.call_args[0][0])
        self.assertFalse(broken.exists())  # cleaned before create

    def test_pip_install_editable_uses_no_deps(self):
        # The env is fully conda-provisioned; pip must only link source, never re-resolve
        # deps from PyPI over the conda packages.
        repo = self.tmp / "repo"
        repo.mkdir()
        with mock.patch.object(env, "_run", return_value=_ok()) as m_run:
            env.pip_install_editable("micromamba", "muagene", repo)
        cmd = m_run.call_args[0][0]
        self.assertIn("--no-deps", cmd)
        self.assertIn("-e", cmd)

    def test_clean_broken_env_prefix_removes_dir_without_conda_meta(self):
        envs_dir = self.tmp / "envs"
        envs_dir.mkdir()
        broken_prefix = envs_dir / "grn"
        broken_prefix.mkdir()
        (broken_prefix / "lib").mkdir()  # has files, no conda-meta

        info_json = '{"envs_dirs": ["' + str(envs_dir) + '"]}'
        with mock.patch.object(env.subprocess, "run",
                               return_value=subprocess.CompletedProcess([], 0, info_json, "")):
            env._clean_broken_env_prefix("micromamba", "grn")

        self.assertFalse(broken_prefix.exists())

    def test_clean_broken_env_prefix_removes_partial_env_with_empty_conda_meta(self):
        # The real-world failure: an interrupted create leaves a prefix with bin/python and
        # an empty conda-meta (no `history`). `conda env list` hides it, so it reads as
        # 'missing' and we try to create — but the directory is there. A conda-meta DIR
        # alone is not proof of health; `conda-meta/history` is. Must be cleaned.
        envs_dir = self.tmp / "envs"
        envs_dir.mkdir()
        broken_prefix = envs_dir / "grn"
        (broken_prefix / "conda-meta").mkdir(parents=True)   # present but EMPTY (no history)
        (broken_prefix / "bin").mkdir()
        (broken_prefix / "bin" / "python").write_text("")

        info_json = '{"envs_dirs": ["' + str(envs_dir) + '"]}'
        with mock.patch.object(env.subprocess, "run",
                               return_value=subprocess.CompletedProcess([], 0, info_json, "")):
            env._clean_broken_env_prefix("micromamba", "grn")

        self.assertFalse(broken_prefix.exists())

    def test_clean_broken_env_prefix_leaves_healthy_env_alone(self):
        envs_dir = self.tmp / "envs"
        envs_dir.mkdir()
        good_prefix = envs_dir / "grn"
        (good_prefix / "conda-meta").mkdir(parents=True)
        (good_prefix / "conda-meta" / "history").write_text("")  # conda's marker of a real env

        info_json = '{"envs_dirs": ["' + str(envs_dir) + '"]}'
        with mock.patch.object(env.subprocess, "run",
                               return_value=subprocess.CompletedProcess([], 0, info_json, "")):
            env._clean_broken_env_prefix("micromamba", "grn")

        self.assertTrue(good_prefix.exists())
        self.assertTrue((good_prefix / "conda-meta" / "history").exists())

    def test_envs_dirs_fallback_when_manager_omits_envs_dirs(self):
        # mamba 2.x's `info --json` has no `envs_dirs` key. Cleanup must still locate the
        # prefix via the derived install-root envs dir, else a failed provision can never
        # self-heal on such a host.
        root = self.tmp / "miniforge3"
        (root / "bin").mkdir(parents=True)
        manager_bin = root / "bin" / "mamba"
        manager_bin.write_text("")
        broken_prefix = root / "envs" / "grn"
        (broken_prefix / "conda-meta").mkdir(parents=True)  # partial: no history

        with mock.patch.object(env.subprocess, "run",
                               return_value=subprocess.CompletedProcess([], 0, "{}", "")), \
             mock.patch.object(env.shutil, "which", return_value=str(manager_bin)):
            dirs = env._conda_envs_dirs("mamba")
            self.assertIn(root / "envs", dirs)
            env._clean_broken_env_prefix("mamba", "grn")

        self.assertFalse(broken_prefix.exists())


class NamedEnvUsableTests(unittest.TestCase):
    """`_named_env_usable` (and thus env_present) must treat a present-but-broken env as
    unusable so it re-provisions, while never false-negativing a healthy or unlocatable env."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.envs = self.tmp / "envs"

    def _prefix(self, name, *, healthy):
        p = self.envs / name
        (p / "conda-meta").mkdir(parents=True)
        if healthy:
            (p / "conda-meta" / "history").write_text("")
        return p

    def test_not_listed_is_unusable(self):
        with mock.patch.object(env, "_conda_env_present", return_value=False):
            self.assertFalse(env._named_env_usable("micromamba", "muagene"))

    def test_listed_and_healthy_is_usable(self):
        self._prefix("muagene", healthy=True)
        with mock.patch.object(env, "_conda_env_present", return_value=True), \
             mock.patch.object(env, "_conda_envs_dirs", return_value=[self.envs]):
            self.assertTrue(env._named_env_usable("micromamba", "muagene"))

    def test_listed_but_broken_prefix_is_unusable(self):
        self._prefix("muagene", healthy=False)  # conda-meta but no history
        with mock.patch.object(env, "_conda_env_present", return_value=True), \
             mock.patch.object(env, "_conda_envs_dirs", return_value=[self.envs]):
            self.assertFalse(env._named_env_usable("micromamba", "muagene"))

    def test_listed_but_unlocatable_prefix_trusts_listing(self):
        # Custom envs_dir we can't see -> don't false-negative a healthy listed env.
        with mock.patch.object(env, "_conda_env_present", return_value=True), \
             mock.patch.object(env, "_conda_envs_dirs", return_value=[self.envs]):
            self.assertTrue(env._named_env_usable("micromamba", "muagene"))


class ReconcileStaleWarningTests(unittest.TestCase):
    """reconcile() emits env_stale_reprovision warning before updating a stale conda env."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        patcher = mock.patch.object(env, "_state_path", return_value=self.tmp / "state.json")
        patcher.start()
        self.addCleanup(patcher.stop)
        self.yaml = self.tmp / "processing.yaml"
        self.yaml.write_text("name: muagene\ndependencies: [scanpy]\n")
        self.lock = self.tmp / "processing.linux-64.lock"
        self.lock.write_text("# platform: linux-64\n@EXPLICIT\nhttps://x/a.conda\n")
        self.sc = SiteConfig(scheduler="slurm", conda_env="grn",
                             environments={"policy": "auto", "cpu": {
                                 "provider": "lock", "definition": str(self.yaml),
                                 "lock": str(self.lock)}})

    def test_stale_warning_emitted_before_reprovision(self):
        with mock.patch.object(env, "_host_conda_subdir", return_value="linux-64"), \
             mock.patch.object(env, "env_status", return_value="stale"), \
             mock.patch.object(env, "provision_env",
                               return_value={"status": "provisioned", "action": "update_from_lock",
                                             "device": "cpu"}), \
             mock.patch.object(env, "validate_env",
                               return_value={"ok": True, "findings": []}):
            out = env.reconcile(self.sc, "/repo", "cpu")
        self.assertTrue(any(f["code"] == "env_stale_reprovision" and f["severity"] == "warning"
                            for f in out["findings"]))

    def test_no_stale_warning_when_env_ok(self):
        with mock.patch.object(env, "_host_conda_subdir", return_value="linux-64"), \
             mock.patch.object(env, "env_status", return_value="ok"), \
             mock.patch.object(env, "validate_env",
                               return_value={"ok": True, "findings": []}):
            out = env.reconcile(self.sc, "/repo", "cpu")
        self.assertFalse(any(f["code"] == "env_stale_reprovision" for f in out["findings"]))


if __name__ == "__main__":
    unittest.main()
