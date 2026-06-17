"""Machine-level infrastructure profile — Execution-MuAgent's record of THIS host.

Top-level rule: Execution-MuAgent owns non-scientific infrastructure. A machine's
*facts* (which env manager, which container runtime, which `module load` provides
singularity, where the GPU image lives, the auto/manual policy, where the
Processing-MuAgent checkout is) are machine-level, not per-run. They used to be
re-typed into every per-run `site.config` (a Processing science deliverable); that
layering smell is fixed here.

`init-machine` writes `~/.muagene/machine.config` once. After that:
  * Execution `provision-env`/`validate-env` can run with NO science site.config —
    they synthesize one from this profile + the committed env manifest (so the
    bootstrap chicken-and-egg disappears: the FIRST command on a fresh machine is an
    Execution command, not a manual `conda env create`).
  * Processing `configure-execution` reads this file to auto-fill machine knobs, so a
    user never re-enters manager/module/image per run.

The two repos cannot import each other, so this file is a stable on-disk YAML
contract. The env-definition *paths* live in exactly one committed place —
``<processing_repo>/workflow/envs/manifest.yaml`` — read by both sides.
"""
from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from . import capabilities
from .monitor import SiteConfig

# Repo root of THIS (Execution-MuAgent) checkout — used by init-machine to install
# the package into the science env so `submit` can spawn `python -m execution_muagent`.
EXECUTION_REPO_ROOT: Path = Path(__file__).resolve().parent.parent

# Conventional default env identities for a fresh machine (overridable via flags).
DEFAULT_CPU_ENV = "muagene"
DEFAULT_GPU_ENV = "muagene-gpu"

# Relative location of the shared env-definition manifest inside Processing-MuAgent.
ENV_MANIFEST_REL = "workflow/envs/manifest.yaml"


def machine_config_path() -> Path:
    return Path(os.path.expanduser("~/.muagene/machine.config"))


@dataclass
class MachineConfig:
    """Per-host infrastructure facts (YAML at ~/.muagene/machine.config)."""
    schema_version: str = "1"
    processing_repo: str | None = None
    manager: str | None = None             # micromamba | mamba | conda
    container_runtime: str | None = None   # apptainer | singularity
    singularity_module: str | None = None  # `module load` name for singularity on HPC
    gpu_image: str | None = None           # machine-local .sif path the image pulls to
    gpu_image_uri: str | None = None       # pinned registry ref the image is PULLED from
    policy: str = "auto"                   # auto | manual (submit-time reconcile policy)
    conda_env: str = DEFAULT_CPU_ENV       # provisioned CPU env name
    gpu_conda_env: str = DEFAULT_GPU_ENV   # provisioned GPU env name (labels the image)
    # Detected at init time, informational (login nodes are not authoritative for GPU).
    scheduler: str | None = None
    gpu_present: bool = False

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "MachineConfig":
        fields = {f for f in cls.__dataclass_fields__}  # type: ignore[attr-defined]
        return cls(**{k: v for k, v in (data or {}).items() if k in fields})


def load_machine_config(path: str | Path | None = None) -> MachineConfig | None:
    """Read ~/.muagene/machine.config. None when absent (fresh, un-bootstrapped host)."""
    import yaml

    p = Path(path) if path else machine_config_path()
    if not p.exists():
        return None
    try:
        data = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    except (OSError, ValueError):
        return None
    return MachineConfig.from_dict(data)


def write_machine_config(cfg: MachineConfig, path: str | Path | None = None) -> Path:
    import yaml

    p = Path(path) if path else machine_config_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(yaml.safe_dump(cfg.to_dict(), default_flow_style=False, sort_keys=False),
                 encoding="utf-8")
    return p


def load_env_manifest(processing_repo: str | Path) -> dict[str, Any]:
    """Load the committed env-definition manifest from the Processing-MuAgent repo.

    This is the single source for the per-device provider + definition/lock/imports
    paths; both agents read it so the path list lives in exactly one file.
    """
    import yaml

    repo = Path(os.path.expanduser(str(processing_repo)))
    mpath = repo / ENV_MANIFEST_REL
    if not mpath.exists():
        raise FileNotFoundError(
            f"env manifest not found at {mpath}. Is --processing-repo a Processing-MuAgent "
            f"checkout? (expected {ENV_MANIFEST_REL} inside it)")
    return yaml.safe_load(mpath.read_text(encoding="utf-8")) or {}


def default_environments_section(processing_repo: str | Path,
                                 cfg: MachineConfig | None = None) -> dict[str, Any]:
    """Build the same `environments:` dict shape Processing writes into site.config,
    sourced from the committed manifest + this machine's profile. Used when there is
    no per-run site.config (init-machine / a bare provision-env)."""
    man = load_env_manifest(processing_repo)
    cfg = cfg or MachineConfig()
    cpu = dict(man.get("cpu") or {})
    gpu = dict(man.get("gpu") or {})
    defaults = dict(man.get("defaults") or {})
    gpu_image = cfg.gpu_image or os.path.expanduser(defaults.get("gpu_image") or "")
    return {
        "manager": cfg.manager,
        "container_runtime": cfg.container_runtime,
        "singularity_module": cfg.singularity_module,
        "policy": cfg.policy or "auto",
        "cpu": {
            "provider": cpu.get("provider") or "lock",
            "definition": cpu.get("definition"),
            "lock": cpu.get("lock"),
            "imports": cpu.get("imports"),
        },
        "gpu": {
            "provider": gpu.get("provider") or "container",
            "definition": gpu.get("definition"),
            "image": gpu_image or None,
            "image_uri": cfg.gpu_image_uri,
            "imports": gpu.get("imports"),
        },
    }


def synthesize_site_config(processing_repo: str | Path, cfg: MachineConfig | None = None,
                           *, device: str = "cpu") -> SiteConfig:
    """A site.config-shaped object built from the machine profile + manifest, so the
    existing resolve_env_spec / provision_env / validate_env path works unchanged
    without a science site.config."""
    cfg = cfg or MachineConfig()
    return SiteConfig(
        scheduler=cfg.scheduler or "local",
        device=device,
        conda_env=cfg.conda_env or DEFAULT_CPU_ENV,
        gpu_conda_env=cfg.gpu_conda_env or DEFAULT_GPU_ENV,
        environments=default_environments_section(processing_repo, cfg),
    )


def detect_machine_config(processing_repo: str | Path, *, manager: str | None = None,
                          container_runtime: str | None = None,
                          singularity_module: str | None = None,
                          gpu_image: str | None = None, gpu_image_uri: str | None = None,
                          policy: str = "auto", conda_env: str | None = None,
                          gpu_conda_env: str | None = None) -> MachineConfig:
    """Probe this host and merge with explicit overrides into a MachineConfig.

    Explicit args win; otherwise capabilities are auto-detected. The GPU image path
    defaults from the manifest's `defaults.gpu_image`.
    """
    caps = capabilities.probe_capabilities()
    man = load_env_manifest(processing_repo)
    default_img = os.path.expanduser((man.get("defaults") or {}).get("gpu_image") or "")
    return MachineConfig(
        processing_repo=str(Path(os.path.expanduser(str(processing_repo))).resolve()),
        manager=manager or caps.get("manager"),
        container_runtime=container_runtime or caps.get("container_runtime"),
        singularity_module=singularity_module,
        gpu_image=gpu_image or default_img or None,
        gpu_image_uri=gpu_image_uri,
        policy=policy or "auto",
        conda_env=conda_env or DEFAULT_CPU_ENV,
        gpu_conda_env=gpu_conda_env or DEFAULT_GPU_ENV,
        scheduler=caps.get("scheduler"),
        gpu_present=bool(caps.get("gpu_present")),
    )
