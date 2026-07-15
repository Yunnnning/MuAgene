"""Local Snakemake invocation — construct + run the snakemake command, and unlock.

`run`/`propose` are the only LOCAL entry points; all cluster execution is owned by
Execution-MuAgent and reached via `submit`. The head-job's own snakemake invocation
(in launch_runner.sh) attaches the cluster profile — this module never does.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import click

from .run_paths import RunPaths

PACKAGE_DIR = Path(__file__).resolve().parent.parent  # Processing-MuAgent/
SNAKEFILE = PACKAGE_DIR / "workflow" / "Snakefile"


def unlock_snakemake(run_dir: Path, config_path: Path) -> None:
    """Run `snakemake --unlock` on the run's workdir (clears a stale lock; no execution)."""
    paths = RunPaths(run_dir)
    env = os.environ.copy()
    env.setdefault("PYTHONPATH", str(PACKAGE_DIR))
    env.setdefault("PMA_REPO_ROOT", str(PACKAGE_DIR))
    paths.snakemake_workdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "-m", "snakemake",
        "-s", str(SNAKEFILE),
        "--directory", str(paths.snakemake_workdir),
        "--unlock",
        "--configfile", str(config_path),
    ]
    click.echo(f"$ {' '.join(cmd)}")
    result = subprocess.run(cmd, env=env, cwd=str(PACKAGE_DIR))
    if result.returncode != 0:
        raise click.ClickException(f"snakemake --unlock exited with {result.returncode}")


def run_snakemake(args: list[str], run_dir: Path) -> None:
    """Invoke snakemake LOCALLY with --cores 1 for reproducibility.

    `run` and `propose` are local-only entry points. All cluster execution is
    owned by Execution-MuAgent and reached via `submit` (which renders + submits
    a supervised head-job) — never through this helper. The head-job's own
    snakemake invocation (in launch_runner.sh) attaches the cluster profile; this
    helper does not.

    Expected args shape from callers: ["--configfile", <path>, <target>].
    """
    env = os.environ.copy()
    env.setdefault("PYTHONPATH", str(PACKAGE_DIR))
    env.setdefault("PMA_REPO_ROOT", str(PACKAGE_DIR))
    paths = RunPaths(run_dir)
    paths.snakemake_workdir.mkdir(parents=True, exist_ok=True)
    env.setdefault("XDG_CACHE_HOME", str(paths.snakemake_workdir / "cache"))
    # Single-thread for reproducibility (UMAP / numba) — unchanged on local;
    # cluster jobs inherit these unless the user overrides in their shell.
    env.setdefault("NUMBA_NUM_THREADS", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("PYTHONHASHSEED", "0")
    if os.environ.get("PMA_AUTO_APPROVE"):
        env["PMA_AUTO_APPROVE"] = os.environ["PMA_AUTO_APPROVE"]

    configfile_path = None
    targets: list[str] = []
    rest: list[str] = []
    it = iter(args)
    for a in it:
        if a == "--configfile":
            configfile_path = next(it, None)
        elif a.startswith("-"):
            rest.append(a)
        else:
            targets.append(a)

    cmd = [
        sys.executable, "-m", "snakemake",
        "-s", str(SNAKEFILE),
        "--directory", str(paths.snakemake_workdir),
        # Rerun only on mtime/missing-output, not params/code/input-set/software-env.
        # This pipeline forces reruns by explicit artifact deletion (executor revise
        # -> _invalidate_qc_downstream), so content/input-set triggers only cause
        # spurious reruns. Mirrors `rerun-triggers: [mtime]` in the cluster profiles.
        "--rerun-triggers", "mtime",
        "--rerun-incomplete", *targets, *rest,
        "--cores", "1",
    ]

    if configfile_path:
        cmd += ["--configfile", configfile_path]
    click.echo(f"$ {' '.join(cmd)}")
    r = subprocess.run(cmd, env=env, cwd=str(PACKAGE_DIR))
    if r.returncode != 0:
        raise click.ClickException(f"snakemake exited with {r.returncode}")
