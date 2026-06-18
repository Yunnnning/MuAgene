"""Environment provisioning + validation — Execution-MuAgent's runtime-infra job.

Processing-MuAgent authors WHAT the science needs: env definitions/locks/`.def`
in its repo, plus the `environments:` section of site.config that says, per device,
the provider (lock | container | yaml) and the in-repo definition. Execution-MuAgent
owns the rest on whatever machine: make it real (`provision_env`), confirm it works
(`validate_env`), record what it installed (a fingerprint marker), and reconcile when
the definition changes (`env_status` → re-provision).

The **fingerprint contract**: the required fingerprint is the provider's identity —
the lock *content* (CPU, conda-lock) or the pinned image *reference* (GPU container,
a registry tag/digest). The provisioned fingerprint is stored in a machine-local
marker. `ok` when they match, `stale` when it changed (a new lock, or a republished
image tag), `missing` when the env/image isn't there. The next submit re-provisions a
stale/missing env automatically (policy=auto).

GPU is **pull-only**: the container image is built + published centrally (out of band,
from `muagene-gpu.def`) and every machine PULLS a pinned image. No machine builds a
container locally — so the fakeroot/subuid build path is gone, and the GPU identity is
the pinned `image_uri`, not the `.def` hash (the `.def` is provenance only).
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import capabilities

# Generous timeouts: an image pull or a from-lock create is a one-time, per-machine,
# per-fingerprint cost. Import validation is quick.
_PULL_TIMEOUT_S = 3600
_CREATE_TIMEOUT_S = 1800
_VALIDATE_TIMEOUT_S = 300


@dataclass
class EnvSpec:
    """Resolved provisioning recipe for one device's env (paths made absolute)."""
    device: str                 # cpu | gpu
    provider: str               # lock | container | yaml
    env_name: str | None        # conda env identity (lock/yaml); also labels the image
    definition: Path | None     # .yaml (lock/yaml) or .def (container build provenance)
    lock: Path | None           # lockfile (lock provider)
    image: Path | None          # local .sif destination (container provider)
    imports: Path | None         # import-check manifest
    image_uri: str | None = None # pinned registry ref to PULL the image from (container)
    singularity_module: str | None = None  # module to `module load` before pull/exec


def _abs(repo_root: Path, value: str | None) -> Path | None:
    if not value:
        return None
    p = Path(os.path.expanduser(value))
    return p if p.is_absolute() else (repo_root / p)


def resolve_env_spec(site_config: Any, repo_root: str | Path, device: str) -> EnvSpec:
    """Build an EnvSpec from site_config.environments + identity, with back-compat
    defaults for site.configs written before the environments contract existed."""
    repo_root = Path(repo_root)
    envs = getattr(site_config, "environments", None) or {}
    section = dict(envs.get(device) or {})
    name = site_config.gpu_conda_env if device == "gpu" else site_config.conda_env
    # Defaults: GPU is a container, CPU is its yaml (until a lock is generated).
    provider = section.get("provider") or ("container" if device == "gpu" else "yaml")
    return EnvSpec(
        device=device,
        provider=provider,
        env_name=name,
        definition=_abs(repo_root, section.get("definition")),
        lock=_abs(repo_root, section.get("lock")),
        image=_abs(repo_root, section.get("image")),
        imports=_abs(repo_root, section.get("imports")),
        image_uri=section.get("image_uri"),
        singularity_module=envs.get("singularity_module"),
    )


# --- fingerprint contract --------------------------------------------------

def _fingerprint_source(spec: EnvSpec) -> Path | None:
    if spec.provider == "lock":
        return spec.lock
    return spec.definition  # yaml provider


def compute_fingerprint(spec: EnvSpec) -> str:
    """The provider's identity — the *required* fingerprint.

    container: the pinned image reference (registry tag/digest), since the machine
    pulls rather than builds — republishing a new tag is what makes the env stale.
    lock/yaml: sha256 of the in-repo source content.
    """
    if spec.provider == "container":
        return f"uri:{spec.image_uri}" if spec.image_uri else ""
    src = _fingerprint_source(spec)
    if not src or not Path(src).exists():
        return ""
    return "sha256:" + hashlib.sha256(Path(src).read_bytes()).hexdigest()


def _state_path() -> Path:
    return Path(os.path.expanduser("~/.muagene/env_state.json"))


def _state_key(spec: EnvSpec) -> str:
    ident = str(spec.image) if spec.provider == "container" else (spec.env_name or "")
    return f"{spec.provider}:{ident}"


def _read_state() -> dict[str, Any]:
    p = _state_path()
    if p.exists():
        try:
            return json.loads(p.read_text())
        except (OSError, ValueError):
            return {}
    return {}


def _record_provisioned(spec: EnvSpec, fingerprint: str) -> None:
    state = _read_state()
    state[_state_key(spec)] = {"fingerprint": fingerprint, "env_name": spec.env_name,
                               "provider": spec.provider}
    p = _state_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")


def _provisioned_fingerprint(spec: EnvSpec) -> str | None:
    entry = _read_state().get(_state_key(spec))
    return entry.get("fingerprint") if entry else None


# --- presence + status -----------------------------------------------------

def _conda_env_present(manager: str, name: str) -> bool:
    try:
        out = subprocess.run([manager, "env", "list"], capture_output=True,
                             text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return False
    for line in out.splitlines():
        tok = line.split()
        if tok and tok[0] == name:
            return True
    return False


def env_present(spec: EnvSpec, manager: str | None) -> bool:
    if spec.provider == "container":
        return bool(spec.image and spec.image.exists())
    return bool(spec.env_name and manager and _conda_env_present(manager, spec.env_name))


def env_status(spec: EnvSpec, manager: str | None) -> str:
    """ok | missing | stale — the reconcile key for the preflight."""
    if not env_present(spec, manager):
        return "missing"
    required = compute_fingerprint(spec)
    provisioned = _provisioned_fingerprint(spec)
    if required and provisioned and required != provisioned:
        return "stale"
    if required and not provisioned:
        # Present but never recorded by us (e.g. a hand-made env). Treat as ok but
        # adopt the fingerprint so future drift is detectable.
        _record_provisioned(spec, required)
    return "ok"


# --- provisioning ----------------------------------------------------------

def provision_env(spec: EnvSpec, site_config: Any, *, manager: str | None = None,
                  container_runtime: str | None = None, force: bool = False) -> dict[str, Any]:
    """Idempotently make the env real for `spec`. No-op when status is ok (unless
    force). Records the fingerprint on success. Returns a structured result.

    CPU = conda-lock create. GPU = PULL a pinned, centrally-published image (never a
    local build). Records the fingerprint on success.
    """
    manager = manager or (getattr(site_config, "environments", {}) or {}).get("manager") \
        or capabilities.detect_manager()
    container_runtime = container_runtime \
        or (getattr(site_config, "environments", {}) or {}).get("container_runtime") \
        or capabilities.detect_container_runtime() or "singularity"

    if not force and env_status(spec, manager) == "ok":
        return {"status": "ok", "action": "noop", "device": spec.device}

    if spec.provider == "container":
        result = _pull_image(spec, container_runtime)
    elif spec.provider == "lock":
        result = _create_from_lock(spec, manager)
    else:
        result = _create_from_yaml(spec, manager)

    if result.get("returncode") == 0:
        _record_provisioned(spec, compute_fingerprint(spec))
        result["status"] = "provisioned"
    else:
        result["status"] = "failed"
    result["device"] = spec.device
    return result


def _pull_image(spec: EnvSpec, runtime: str) -> dict[str, Any]:
    """Pull the GPU container from its pinned registry reference. No machine ever
    builds the image locally (the build + push is a central, out-of-band maintainer
    step from muagene-gpu.def) — so there is no `--fakeroot`/subuid requirement here."""
    if not spec.image or not spec.image_uri:
        return {"returncode": 1, "action": "pull_image", "code": "gpu_image_unavailable",
                "stderr": "GPU image needs a local `image` path and a pinned `image_uri` to pull "
                          "from (e.g. docker://<registry>/muagene-gpu:<tag>). No machine builds "
                          "the container locally; set environments.gpu.image_uri (configure-execution "
                          "--gpu-image-uri / init-machine --gpu-image-uri)."}
    spec.image.parent.mkdir(parents=True, exist_ok=True)
    mod = f"module load {spec.singularity_module}; " if spec.singularity_module else ""
    # --force overwrites a stale .sif when the published tag/digest changed.
    cmd = (f"{mod}{runtime} pull --force "
           f"{shlex.quote(str(spec.image))} {shlex.quote(str(spec.image_uri))}")
    proc = _bash_login(cmd, _PULL_TIMEOUT_S)
    return {"returncode": proc.returncode, "action": "pull_image", "image": str(spec.image),
            "image_uri": spec.image_uri, "stdout": proc.stdout[-2000:], "stderr": proc.stderr[-2000:]}


def _clean_broken_env_prefix(manager: str, env_name: str) -> None:
    """Remove a non-conda directory at the env prefix so conda can create there.

    This happens when a previous provision attempt removed conda-meta but failed before
    completing the recreation (e.g. an NFS busy file blocked the removal). Without
    cleanup, conda create raises 'Non-conda folder exists at prefix'. Silently skips if
    the prefix is healthy, cannot be located, or the rmtree fails (conda's own error
    will then explain what happened).
    """
    try:
        out = subprocess.run([manager, "info", "--json"], capture_output=True,
                             text=True, timeout=30)
        info = json.loads(out.stdout)
        for d in (info.get("envs_dirs") or []):
            prefix = Path(d) / env_name
            if prefix.exists() and not (prefix / "conda-meta").exists():
                import shutil
                shutil.rmtree(prefix, ignore_errors=True)
                return
    except Exception:
        pass


def _create_from_lock(spec: EnvSpec, manager: str | None) -> dict[str, Any]:
    if not manager or not spec.env_name or not spec.lock:
        return {"returncode": 1, "action": "create_from_lock",
                "stderr": "lock provider needs a manager, env_name, and lock file"}
    if _conda_env_present(manager, spec.env_name):
        # Stale env: update packages in place — no deletion, so a failed update leaves
        # the previous working env intact (no broken-env-on-NFS-failure risk).
        cmd = [manager, "install", "-y", "-n", spec.env_name, "--file", str(spec.lock)]
        action = "update_from_lock"
    else:
        # Missing env: clean up any broken prefix left by a prior failed provision,
        # then create fresh.
        _clean_broken_env_prefix(manager, spec.env_name)
        cmd = [manager, "create", "-y", "-n", spec.env_name, "--file", str(spec.lock)]
        action = "create_from_lock"
    proc = _run(cmd, _CREATE_TIMEOUT_S)
    return {"returncode": proc.returncode, "action": action,
            "env": spec.env_name, "stdout": proc.stdout[-2000:], "stderr": proc.stderr[-2000:]}


def _create_from_yaml(spec: EnvSpec, manager: str | None) -> dict[str, Any]:
    if not manager or not spec.definition:
        return {"returncode": 1, "action": "create_from_yaml",
                "stderr": "yaml provider needs a manager and a definition yaml"}
    proc = _run([manager, "env", "create", "-f", str(spec.definition)], _CREATE_TIMEOUT_S)
    return {"returncode": proc.returncode, "action": "create_from_yaml",
            "env": spec.env_name, "stdout": proc.stdout[-2000:], "stderr": proc.stderr[-2000:]}


def pip_install_editable(manager: str | None, env_name: str | None,
                         repo: str | Path) -> dict[str, Any]:
    """`<manager> run -n <env> pip install -e <repo>` — used by init-machine to layer
    the agent packages onto the provisioned CPU env. Returns a structured result."""
    repo_path = Path(os.path.expanduser(str(repo)))
    if not manager or not env_name:
        return {"returncode": 1, "action": "pip_install",
                "stderr": "pip install needs a manager and an env_name", "repo": str(repo_path)}
    if not repo_path.exists():
        return {"returncode": 1, "action": "pip_install",
                "stderr": f"repo not found: {repo_path}", "repo": str(repo_path)}
    proc = _run([manager, "run", "-n", env_name, "pip", "install", "-e", str(repo_path)],
                _CREATE_TIMEOUT_S)
    return {"returncode": proc.returncode, "action": "pip_install", "repo": str(repo_path),
            "stdout": proc.stdout[-2000:], "stderr": proc.stderr[-2000:]}


# --- validation ------------------------------------------------------------

def _host_conda_subdir() -> str:
    """This host's conda platform subdir (linux-64, osx-arm64, …)."""
    sysname = platform.system().lower()
    machine = platform.machine().lower()
    if sysname == "linux":
        return "linux-aarch64" if machine in ("aarch64", "arm64") else "linux-64"
    if sysname == "darwin":
        return "osx-arm64" if machine in ("arm64", "aarch64") else "osx-64"
    if sysname == "windows":
        return "win-64"
    return f"{sysname}-{machine}"


def _lock_header(lock: Path, key: str) -> str | None:
    """Read a `# <key>: <value>` marker from a lockfile's leading comment lines."""
    try:
        for line in lock.read_text().splitlines()[:12]:
            low = line.lower()
            if low.startswith(f"# {key}:"):
                return line.split(":", 1)[1].strip()
    except OSError:
        return None
    return None


def _lock_platform(lock: Path) -> str | None:
    """The platform a conda-lock lockfile was generated for (its `# platform:` header)."""
    return _lock_header(lock, "platform")


def _file_sha256(path: Path) -> str | None:
    try:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    except OSError:
        return None


def _lock_preflight_findings(spec: EnvSpec) -> list[dict[str, str]]:
    """Fail-loud guards for the lock (CPU) provider — run before any create/validate:

    * platform_unsupported: the lock is for a different OS/arch than this host. MuAgene's
      CPU env is linux-only; never silently solve the YAML on the wrong platform.
    * lock_stale_vs_yaml: the lock's recorded `# source-sha256:` (the hash of the YAML it
      was generated from, written by `regenerate-locks`) no longer matches the current
      processing.yaml — the lock is out of date and must be regenerated. Content-hash,
      not mtime: git does not preserve mtimes, so an mtime check would false-fire on a
      fresh clone; the hash is stable across clones.
    """
    out: list[dict[str, str]] = []
    if spec.provider != "lock" or not spec.lock:
        return out
    lock = Path(spec.lock)
    lock_subdir = _lock_platform(lock) if lock.exists() else None
    host = _host_conda_subdir()
    if lock_subdir and lock_subdir != host:
        out.append({"severity": "error", "code": "platform_unsupported",
            "message": f"CPU env lock is {lock_subdir} but this host is {host}. MuAgene's CPU env is "
                       f"linux-only — use a {lock_subdir} host (or run via a container), never a "
                       f"silent cross-platform solve."})
    if spec.definition and Path(spec.definition).exists() and lock.exists():
        recorded = _lock_header(lock, "source-sha256")
        current = _file_sha256(Path(spec.definition))
        # Only flag when the lock RECORDS a source hash and it mismatches. A lock with no
        # marker predates the convention (can't verify) — don't false-fire; the next
        # `regenerate-locks` stamps it.
        if recorded and current and recorded != current:
            out.append({"severity": "error", "code": "lock_stale_vs_yaml",
                "message": f"{Path(spec.definition).name} changed since {lock.name} was generated "
                           f"(source hash mismatch); the YAML is the source of truth, the lock is what "
                           f"installs. Regenerate it: `Processing-MuAgent regenerate-locks`, then commit."})
    return out


def _read_imports(spec: EnvSpec) -> list[str]:
    if not spec.imports or not Path(spec.imports).exists():
        return []
    mods: list[str] = []
    for line in Path(spec.imports).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            mods.append(line)
    return mods


def validate_env(spec: EnvSpec, site_config: Any, *, manager: str | None = None,
                 container_runtime: str | None = None) -> dict[str, Any]:
    """Confirm the env is present and imports its declared modules. Returns
    {ok: bool, findings: [...]}. A missing module is an error (never silently
    degrade); a CUDA-unavailable import on a non-GPU host is a warning (the real
    job runs on the GPU node)."""
    manager = manager or (getattr(site_config, "environments", {}) or {}).get("manager") \
        or capabilities.detect_manager()
    container_runtime = container_runtime \
        or (getattr(site_config, "environments", {}) or {}).get("container_runtime") \
        or capabilities.detect_container_runtime() or "singularity"

    # Lock-provider guards first (platform + freshness) — a wrong-platform or stale lock
    # is a problem regardless of whether the env happens to exist.
    findings: list[dict[str, str]] = _lock_preflight_findings(spec)
    if any(f["code"] == "platform_unsupported" for f in findings):
        return {"ok": False, "findings": findings}
    if not env_present(spec, manager):
        what = str(spec.image) if spec.provider == "container" else (spec.env_name or "<unnamed>")
        findings.append({"severity": "error", "code": "env_missing",
                "message": f"{spec.device} env not provisioned: {what}. Run `provision-env --device {spec.device}`."})
        return {"ok": False, "findings": findings}

    mods = _read_imports(spec)
    if mods:
        code = "import " + ", ".join(mods)
        proc = _import_check(spec, code, manager, container_runtime)
        if proc.returncode != 0:
            err = (proc.stderr or proc.stdout or "")
            cuda_unavail = any(s in err for s in ("CUDARuntimeError", "no CUDA-capable device",
                                                  "cudaErrorNoDevice", "libcuda", "CUDADriver"))
            if cuda_unavail and spec.device == "gpu" and not capabilities.gpu_present():
                findings.append({"severity": "warning", "code": "gpu_import_needs_node",
                    "message": f"{spec.device} env present; GPU imports could not be fully verified "
                               "on this non-GPU host (they run on the GPU node at job time)."})
            else:
                findings.append({"severity": "error", "code": "import_failed",
                    "message": f"{spec.device} env failed to import {mods}: {err.strip()[-400:]}"})

    fp_req, fp_prov = compute_fingerprint(spec), _provisioned_fingerprint(spec)
    if fp_req and fp_prov and fp_req != fp_prov:
        findings.append({"severity": "warning", "code": "env_stale",
            "message": f"{spec.device} env fingerprint drifted (definition changed); re-provision to refresh."})

    ok = not any(f["severity"] == "error" for f in findings)
    return {"ok": ok, "findings": findings}


def _import_check(spec: EnvSpec, code: str, manager: str | None, runtime: str):
    if spec.provider == "container":
        mod = f"module load {spec.singularity_module}; " if spec.singularity_module else ""
        cmd = f"{mod}{runtime} exec {shlex.quote(str(spec.image))} python -c {shlex.quote(code)}"
        return _bash_login(cmd, _VALIDATE_TIMEOUT_S)
    return _run([manager, "run", "-n", spec.env_name, "python", "-c", code], _VALIDATE_TIMEOUT_S)


# --- high-level reconcile (used by the execute-spec preflight) --------------

def reconcile(site_config: Any, repo_root: str | Path, device: str) -> dict[str, Any]:
    """Preflight: ensure the device's env is provisioned + valid before a submit.

    policy=auto re-provisions a missing/stale env; policy=manual fails loud with the
    exact command. Returns {ok, status, provision, validation, findings}.
    """
    spec = resolve_env_spec(site_config, repo_root, device)
    manager = (getattr(site_config, "environments", {}) or {}).get("manager") \
        or capabilities.detect_manager()
    policy = (getattr(site_config, "environments", {}) or {}).get("policy", "auto")

    # Fail loud before creating an env from a wrong-platform or stale lock — never
    # auto-provision over a lock that no longer matches its source YAML / this host.
    pre = _lock_preflight_findings(spec)
    if any(f["severity"] == "error" for f in pre):
        return {"device": device, "status": "blocked", "policy": policy,
                "provision": None, "validation": None, "findings": pre, "ok": False}

    status = env_status(spec, manager)
    out: dict[str, Any] = {"device": device, "status": status, "policy": policy,
                           "provision": None, "validation": None, "findings": []}
    if status == "stale":
        out["findings"].append({"severity": "warning", "code": "env_stale_reprovision",
            "message": (f"{device} env lock fingerprint changed (e.g. git branch switch "
                        f"or lock update); updating env '{spec.env_name}' in place — "
                        f"may take a few minutes.")})
    if status != "ok":
        if policy == "manual":
            out["ok"] = False
            out["findings"] = [{"severity": "error", "code": f"env_{status}",
                "message": f"{device} env is {status}; run "
                           f"`Execution-MuAgent provision-env --device {device} --site-config <path>` first."}]
            return out
        out["provision"] = provision_env(spec, site_config, manager=manager)
        if out["provision"].get("status") == "failed":
            out["ok"] = False
            out["findings"] = [{"severity": "error", "code": "provision_failed",
                "message": f"{device} env provisioning failed: "
                           f"{out['provision'].get('stderr', '')[-400:]}"}]
            return out

    validation = validate_env(spec, site_config, manager=manager)
    out["validation"] = validation
    out["findings"] = out["findings"] + validation["findings"]
    out["ok"] = validation["ok"]
    return out


# --- subprocess helpers ----------------------------------------------------

def _run(cmd: list[str], timeout_s: int):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
    except (OSError, subprocess.SubprocessError) as exc:
        return subprocess.CompletedProcess(cmd, 1, "", str(exc))


def _bash_login(cmd: str, timeout_s: int):
    # `bash -lc` so HPC `module` is available for `module load singularityce`.
    try:
        return subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True, timeout=timeout_s)
    except (OSError, subprocess.SubprocessError) as exc:
        return subprocess.CompletedProcess(cmd, 1, "", str(exc))
