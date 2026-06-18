"""GPU/device fields in site.config + hpc.env, and SLURM GPU capability parsing."""
from pathlib import Path
import tempfile
import unittest

from executor import hpc


class ParseSlurmGpuTests(unittest.TestCase):
    def test_detects_gpu_partitions_and_suggests_single_gpu_gres(self):
        # Mirrors real `sinfo -h -o '%P|%G'` from the lab cluster.
        text = (
            "nice*|(null)\n"
            "nice*|gpu:A5000:4(S:0-1)\n"
            "cpu|(null)\n"
            "gpu|gpu:A5000:4(S:0-1)\n"
        )
        parts, gres = hpc._parse_slurm_gpu(text)
        self.assertEqual(set(parts), {"nice", "gpu"})
        self.assertEqual(gres, "gpu:A5000:1")  # 4 -> 1, type preserved

    def test_generic_gres_when_type_unnamed(self):
        parts, gres = hpc._parse_slurm_gpu("gpu|gpu:4\n")
        self.assertEqual(parts, ["gpu"])
        self.assertEqual(gres, "gpu:1")

    def test_none_when_no_gpu(self):
        parts, gres = hpc._parse_slurm_gpu("cpu|(null)\nhmem|(null)\n")
        self.assertEqual(parts, [])
        self.assertIsNone(gres)


class SiteConfigGpuRoundTripTests(unittest.TestCase):
    _SLURM_GPU = {
        "slurm_partition": "cpu", "slurm_account": "vaquerizas",
        "conda_env": "muagene", "gpu_conda_env": "muagene-gpu",
        "device": "gpu", "slurm_gpu_partition": "gpu",
        "slurm_gpu_gres": "gpu:A5000:1", "resources_scale": "1",
    }

    def test_write_site_config_includes_gpu_and_device(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = Path(tmp) / "site.config"
            hpc.write_site_config(site, mode="slurm", settings=dict(self._SLURM_GPU))
            cfg = hpc.load_site_config(site)
            self.assertEqual(cfg["device"], "gpu")
            self.assertEqual(cfg["slurm"]["gpu_partition"], "gpu")
            self.assertEqual(cfg["slurm"]["gpu_gres"], "gpu:A5000:1")
            self.assertEqual(cfg["common"]["gpu_conda_env"], "muagene-gpu")

    def test_hpc_env_exports_gpu_vars(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = Path(tmp) / "site.config"
            hpc.write_site_config(site, mode="slurm", settings=dict(self._SLURM_GPU))
            envsh = Path(tmp) / "hpc.env"
            hpc.write_hpc_env(envsh, site)
            text = envsh.read_text()
            self.assertIn("export PMA_DEVICE='gpu'", text)
            self.assertIn("export PMA_SLURM_GPU_PARTITION='gpu'", text)
            self.assertIn("export PMA_SLURM_GPU_GRES='gpu:A5000:1'", text)
            self.assertIn("export PMA_CONDA_ENV_GPU='muagene-gpu'", text)
            # run_dir is a runtime path derived live by launch_runner.sh — it must NOT be
            # projected into hpc.env (a static site.config snapshot would drift).
            self.assertNotIn("PMA_RUN_DIR", text)

    def test_scratch_round_trips_to_pma_gpu_bind(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = Path(tmp) / "site.config"
            settings = dict(self._SLURM_GPU, scratch="/scratch/fast")
            hpc.write_site_config(site, mode="slurm", settings=settings)
            cfg = hpc.load_site_config(site)
            self.assertEqual(cfg["common"]["scratch"], "/scratch/fast")
            envsh = Path(tmp) / "hpc.env"
            hpc.write_hpc_env(envsh, site)
            self.assertIn("export PMA_GPU_BIND='/scratch/fast'", envsh.read_text())

    def test_no_scratch_omits_pma_gpu_bind(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = Path(tmp) / "site.config"
            hpc.write_site_config(site, mode="slurm", settings=dict(self._SLURM_GPU))
            envsh = Path(tmp) / "hpc.env"
            hpc.write_hpc_env(envsh, site)
            self.assertNotIn("PMA_GPU_BIND", envsh.read_text())

    def test_pbs_gpu_select_extra_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            site = Path(tmp) / "site.config"
            hpc.write_site_config(site, mode="pbs", settings={
                "pbs_queue": "workq", "pbs_project": "vaquerizas",
                "conda_env": "muagene", "gpu_conda_env": "muagene-gpu",
                "device": "gpu", "pbs_gpu_select_extra": "ngpus=1",
                "pbs_gpu_queue": "gpuq", "resources_scale": "1",
            })
            cfg = hpc.load_site_config(site)
            self.assertEqual(cfg["pbs"]["gpu_select_extra"], "ngpus=1")
            self.assertEqual(cfg["pbs"]["gpu_queue"], "gpuq")


class ConfigureExecutionLocalGpuTests(unittest.TestCase):
    def test_local_mode_rejects_device_gpu(self):
        import yaml
        from click.testing import CliRunner
        from executor import cli

        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "run"
            draft = Path(tmp) / "draft.yaml"
            draft.write_text(yaml.safe_dump({"run_dir": str(run_dir)}))
            runner = CliRunner()
            assert runner.invoke(cli.main, ["init", "--config", str(draft)]).exit_code == 0
            cfg = run_dir / "deliverables" / "plan" / "config" / "run.yaml"
            res = runner.invoke(cli.main, [
                "configure-execution", "--config", str(cfg),
                "--mode", "local", "--device", "gpu", "--confirmed-by-user",
            ])
            self.assertNotEqual(res.exit_code, 0)
            self.assertIn("cluster-only", res.output)


class ConfigureExecutionGpuImageUriTests(unittest.TestCase):
    """SLURM --device gpu must fail loud at configure time when no image_uri is
    resolvable — container is the only SLURM GPU provider and the image is PULLED from
    the pinned reference, so a missing one is caught now, not at provision/submit."""

    def _init_run(self, tmp):
        import yaml
        from click.testing import CliRunner
        from executor import cli
        run_dir = Path(tmp) / "run"
        draft = Path(tmp) / "draft.yaml"
        draft.write_text(yaml.safe_dump({"run_dir": str(run_dir)}))
        runner = CliRunner()
        assert runner.invoke(cli.main, ["init", "--config", str(draft)]).exit_code == 0
        cfg = run_dir / "deliverables" / "plan" / "config" / "run.yaml"
        return runner, cli, cfg

    def test_slurm_gpu_without_image_uri_fails_loud(self):
        from unittest import mock
        with tempfile.TemporaryDirectory() as tmp:
            runner, cli, cfg = self._init_run(tmp)
            # No image_uri anywhere: machine.config empty + env var blanked.
            with mock.patch("executor.hpc.load_machine_config", return_value={}), \
                 mock.patch.dict("os.environ", {"PMA_GPU_IMAGE_URI": ""}):
                res = runner.invoke(cli.main, [
                    "configure-execution", "--config", str(cfg),
                    "--mode", "slurm", "--slurm-partition", "cpu",
                    "--slurm-account", "vaquerizas", "--conda-env", "muagene",
                    "--device", "gpu", "--gpu-gres", "gpu:A5000:1", "--confirmed-by-user",
                ])
            self.assertNotEqual(res.exit_code, 0)
            self.assertIn("--gpu-image-uri", res.output)

    def test_slurm_gpu_with_image_uri_succeeds(self):
        from unittest import mock
        with tempfile.TemporaryDirectory() as tmp:
            runner, cli, cfg = self._init_run(tmp)
            with mock.patch("executor.hpc.load_machine_config", return_value={}), \
                 mock.patch.dict("os.environ", {"PMA_GPU_IMAGE_URI": ""}):
                res = runner.invoke(cli.main, [
                    "configure-execution", "--config", str(cfg),
                    "--mode", "slurm", "--slurm-partition", "cpu",
                    "--slurm-account", "vaquerizas", "--conda-env", "muagene",
                    "--device", "gpu", "--gpu-gres", "gpu:A5000:1",
                    "--gpu-image-uri", "docker://example/muagene-gpu:test",
                    "--confirmed-by-user",
                ])
            self.assertEqual(res.exit_code, 0, res.output)


if __name__ == "__main__":
    unittest.main()
