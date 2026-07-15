"""Golden test: the LOCAL Snakemake invocation must stay byte-for-byte stable.

`run`/`propose` go through snakemake_runner.run_snakemake. A drift in the constructed
argv or env (cores, rerun-triggers, configfile, reproducibility env) would silently change
local pipeline behaviour, so we pin both here by capturing the subprocess call.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from unittest import mock

from executor import snakemake_runner as sr
from executor.run_paths import RunPaths


def _capture(args: list[str]):
    with tempfile.TemporaryDirectory() as tmp:
        run_dir = Path(tmp)
        captured = {}

        def fake_run(cmd, env=None, cwd=None):
            captured["cmd"] = cmd
            captured["env"] = env
            captured["cwd"] = cwd
            return mock.Mock(returncode=0)

        with mock.patch.object(sr.subprocess, "run", side_effect=fake_run):
            sr.run_snakemake(args, run_dir)
        captured["workdir"] = str(RunPaths(run_dir).snakemake_workdir)
        return captured


def test_run_snakemake_builds_expected_argv():
    cap = _capture(["--configfile", "/runs/x/run.yaml", "all"])
    cmd = cap["cmd"]
    assert cmd == [
        sys.executable, "-m", "snakemake",
        "-s", str(sr.SNAKEFILE),
        "--directory", cap["workdir"],
        "--rerun-triggers", "mtime",
        "--rerun-incomplete", "all",
        "--cores", "1",
        "--configfile", "/runs/x/run.yaml",
    ]
    assert cap["cwd"] == str(sr.PACKAGE_DIR)


def test_run_snakemake_sets_reproducibility_env():
    cap = _capture(["--configfile", "/runs/x/run.yaml", "all"])
    env = cap["env"]
    assert env["NUMBA_NUM_THREADS"] == "1"
    assert env["OMP_NUM_THREADS"] == "1"
    assert env["PYTHONHASHSEED"] == "0"
    assert env["PYTHONPATH"].startswith(str(sr.PACKAGE_DIR)) or env["PYTHONPATH"] == str(sr.PACKAGE_DIR)


def test_run_snakemake_separates_targets_flags_configfile():
    # extra flags pass through; multiple targets preserved; configfile pulled to the end.
    cap = _capture(["--configfile", "/c.yaml", "s3_doublets_execute", "--dry-run"])
    cmd = cap["cmd"]
    assert "s3_doublets_execute" in cmd
    assert "--dry-run" in cmd
    assert cmd[-2:] == ["--configfile", "/c.yaml"]
