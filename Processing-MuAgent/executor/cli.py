"""Processing-MuAgent CLI — thin command-entry layer over the executor modules.

Each command resolves the run dir, then delegates to a focused module:
  - pipeline.py          — stage topology (order, per-branch membership, aliases)
  - cleanup.py           — cleanup/reset policy (QC + S4–S8 intermediates)
  - revision.py          — `revise` + the QC-invalidation cascade
  - gates.py             — context / execution-mode / marker-gene gates + approval seeding
  - snakemake_runner.py  — local Snakemake invocation + unlock
  - reporting.py         — `status` / `hpc-status` rendering
The moved helpers are re-bound below under their original (underscore) names so the
public/internal API and test patch targets (e.g. cli._snakemake) are unchanged.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import click
import yaml

from . import approval, context as _ctx, hpc, plan_review as _pr, provenance, reporting, specs as _specs, stage_progress as _sp
from .log import log_event
from .run_paths import RunPaths
from .pipeline import (
    HUMAN_CHECKPOINTS as HUMAN_CHECKPOINT_STAGES,
    canonical_stage as _canonical_stage,
    display_stage as _display_stage,
)
from .cleanup import (
    cleanup_qc_intermediates as _cleanup_qc_intermediates,
    cleanup_process_intermediates as _cleanup_process_intermediates,
    s8_outputs_valid as _s8_outputs_valid,
    ensure_process_markers as _ensure_process_markers,
)
from .revision import (
    qc_downstream_targets as _qc_downstream_targets,
    invalidate_qc_downstream as _invalidate_qc_downstream,
    revise_dry_run as _revise_dry_run,
    regenerate_plan_deliverables as _regenerate_plan_deliverables,
)
from .gates import (
    apply_marker_gene_ack as _apply_marker_gene_ack,
    resolve_marker_gene_gate as _resolve_marker_gene_gate,
    enforce_context_gate as _enforce_context_gate,
    enforce_execution_mode_gate as _enforce_execution_mode_gate,
    infer_submit_target as _infer_submit_target,
    missing_approvals as _missing_approvals,
    seed_approvals as _seed_approvals,
    prepare_submit_approvals as _prepare_submit_approvals,
)
from .snakemake_runner import (
    run_snakemake as _snakemake,
    unlock_snakemake as _unlock_snakemake,
)


EXECUTOR_CHOICE = click.Choice(["local", "slurm"])
# Cluster-only executor for `submit` — `run` is local-only, so `local` is not a
# valid submit target (all cluster execution is owned by Execution-MuAgent).
CLUSTER_EXECUTOR_CHOICE = click.Choice(["slurm"])


def _resolve_run_dir(config_path: Path | str) -> Path:
    with Path(config_path).open() as f:
        cfg = yaml.safe_load(f) or {}
    rd = cfg.get("run_dir")
    if not rd:
        raise click.ClickException("run.yaml must set 'run_dir'")
    return Path(rd).expanduser().resolve()


@click.group()
def main() -> None:
    """Processing-MuAgent: multiome preprocessing subagent (stops after per-modality UMAP)."""


@main.command()
@click.option("--config", "config_path", required=True, type=click.Path())
def init(config_path: str) -> None:
    """Initialize a run directory.

    Creates the `internal/` and `deliverables/` scaffolds, copies the user's
    config into its canonical user-facing location `deliverables/plan/config/run.yaml`,
    and writes the Biological Context Report template into
    `deliverables/plan/config/biological_context.md`.
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    paths.ensure()
    # Config goes to its canonical deliverable location (Snakemake will read it
    # from there via --configfile — no separate internal copy).
    shutil.copy(config_path, paths.run_yaml)
    _ctx.write_template(paths.biological_context_md)
    click.echo(f"Initialized {run_dir}")
    click.echo(f"Fill {paths.biological_context_md} (optional but recommended).")


@main.command()
@click.argument("stage")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def propose(stage: str, config_path: str) -> None:
    """Run the <stage>_propose rule (local — propose rules are localrules)."""
    stage = _canonical_stage(stage)
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    _snakemake(["--configfile", str(paths.run_yaml), f"{stage}_propose"], run_dir)


@main.command()
@click.argument("stage")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--note", default="")
@click.option("--defer-marker-genes", "defer_marker_genes", is_flag=True,
              help="plan_review only: record an explicit choice to check marker "
                   "genes at QC review instead of now.")
@click.option("--skip-marker-genes", "skip_marker_genes", is_flag=True,
              help="plan_review only: record an explicit choice to decline the "
                   "before/after-ambient marker gene expression check.")
def approve(stage: str, config_path: str, note: str,
            defer_marker_genes: bool, skip_marker_genes: bool) -> None:
    """Write internal/checkpoints/<stage>.approved to unblock <stage>_execute."""
    stage = _canonical_stage(stage)
    run_dir = _resolve_run_dir(config_path)
    if stage == "plan_review":
        _resolve_marker_gene_gate(
            run_dir, defer=defer_marker_genes, skip=skip_marker_genes)
    elif defer_marker_genes or skip_marker_genes:
        raise click.ClickException(
            "--defer-marker-genes / --skip-marker-genes apply only to "
            "`approve plan_review`.")
    approval.approve(run_dir, stage, note=note)
    log_event(run_dir, {"stage": stage, "event": "approved", "note": note})
    if stage == "post_qc_review":
        deleted = _cleanup_qc_intermediates(run_dir)
        if deleted:
            log_event(run_dir, {"stage": "post_qc_review", "event": "qc_cleanup",
                                 "deleted": deleted})
            click.echo(f"Cleaned up {len(deleted)} intermediate QC object(s).")
    click.echo(f"Approved {_display_stage(stage)}")


@main.command(name="finish-cleanup")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def finish_cleanup(config_path: str) -> None:
    """Delete the large S4–S8 intermediate working files after the run is complete.

    Run this once S8 has produced the final processed deliverable. It first VALIDATES
    that the S8 output exists and is non-empty; if not, it refuses and leaves every
    intermediate in place so the pipeline can still resume from an intermediate stage.
    On success it backfills any missing durable stage markers (so `status` keeps
    reporting S4–S8 done) and removes the working h5ads/sidecars — content-duplicates
    of the processed deliverable. The deletions are not declared Snakemake outputs, so
    a later `submit --target all` does not re-run S4–S8.
    """
    run_dir = _resolve_run_dir(config_path)
    ok, problems = _s8_outputs_valid(run_dir)
    if not ok:
        raise click.ClickException(
            "S8 output not present/valid — refusing finish-cleanup so the run can "
            "resume from an intermediate step:\n  " + "\n  ".join(problems))
    backfilled = _ensure_process_markers(run_dir)
    if backfilled:
        log_event(run_dir, {"stage": "finish_cleanup", "event": "markers_backfilled",
                            "written": backfilled})
    deleted = _cleanup_process_intermediates(run_dir)
    if deleted:
        log_event(run_dir, {"stage": "finish_cleanup", "event": "process_cleanup",
                            "deleted": deleted})
        click.echo(f"Cleaned up {len(deleted)} S4–S8 intermediate object(s).")
    else:
        click.echo("No S4–S8 intermediates to clean.")


@main.command(name="qc-cleanup")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def qc_cleanup(config_path: str) -> None:
    """Delete the large QC/ingest intermediate working files of an approved run.

    This is the same cleanup `approve post_qc_review` runs automatically; expose it
    standalone so disk can be reclaimed on a run that was approved earlier (e.g. to
    apply an expanded cleanup set retroactively). REQUIRES `post_qc_review` to be
    approved — refuses otherwise, since the deleted caches (rna_ingest.h5ad,
    rna_decontaminated.h5ad, the QC matrices) are still needed while QC is in review.
    The durable markers survive, so nothing re-runs; deliverables are untouched.
    """
    run_dir = _resolve_run_dir(config_path)
    if not approval.is_approved(run_dir, "post_qc_review"):
        raise click.ClickException(
            "post_qc_review is not approved — refusing qc-cleanup. These caches are "
            "needed while QC is in review (and the before/after-ambient marker-gene "
            "check reads rna_decontaminated.h5ad). Approve QC first.")
    deleted = _cleanup_qc_intermediates(run_dir)
    if deleted:
        log_event(run_dir, {"stage": "post_qc_review", "event": "qc_cleanup",
                            "deleted": deleted})
        click.echo(f"Cleaned up {len(deleted)} intermediate QC object(s).")
    else:
        click.echo("No QC intermediates to clean.")


@main.command(name="declare-branch")
@click.argument("branch", type=click.Choice(["paired", "unpaired", "rna_only", "atac_only"]))
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def declare_branch(branch: str, config_path: str) -> None:
    """Declare the workflow branch up front (user assertion).

    Writes `plan.workflow_branch_declared` to parameters.yaml with source=user.
    S0 will confirm this matches its own detection, or raise with a clear diff.
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    paths.ensure()
    provenance.set_param(
        str(paths.parameters_yaml),
        "plan.workflow_branch_declared", branch,
        source="user", confidence="high",
        rationale=f"Declared via `executor declare-branch {branch}`.",
    )
    log_event(run_dir, {"stage": "declare_branch", "event": "declared", "branch": branch})
    click.echo(f"Declared workflow_branch={branch!r}; S0 will confirm at ingest time.")


@main.command(name="hpc-info")
def hpc_info() -> None:
    """Probe the login node for scheduler queues/partitions and current PMA_* env."""
    import json
    info = hpc.discover_site()
    click.echo(json.dumps(info, indent=2, sort_keys=True))


@main.command(name="configure-execution")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--mode", "mode", required=True, type=EXECUTOR_CHOICE,
              help="Execution backend: local or slurm.")
@click.option("--confirmed-by-user/--not-confirmed", "confirmed_by_user",
              default=False,
              help="Record that the USER explicitly confirmed this execution mode "
                   "(local vs HPC). `run`/`submit` refuse to launch any compute job "
                   "until this is set. Never pass it on the user's behalf without "
                   "having actually confirmed.")
@click.option("--slurm-partition", default=None, help="SLURM partition (PMA_SLURM_PARTITION).")
@click.option("--slurm-account", default=None, help="SLURM account (PMA_SLURM_ACCOUNT).")
@click.option("--resources-scale", default=None, type=float,
              help="Memory/walltime scale factor (PMA_RESOURCES_SCALE).")
@click.option("--conda-env", default=None, help="Conda env name for cluster jobs (PMA_CONDA_ENV).")
@click.option("--device", type=click.Choice(["cpu", "gpu"]), default="cpu", show_default=True,
              help="Compute device for GPU-capable stages (compute.device). 'gpu' is cluster-only "
                   "(--mode slurm + submit): routes those stages to the GPU partition/gres and "
                   "container; local mode (--mode local) is always CPU.")
@click.option("--gpu-partition", default=None,
              help="SLURM GPU partition for GPU-capable stages (PMA_SLURM_GPU_PARTITION).")
@click.option("--gpu-gres", default=None,
              help="SLURM GPU gres request, e.g. gpu:A5000:1 (PMA_SLURM_GPU_GRES). Run hpc-info to discover.")
@click.option("--gpu-conda-env", default=None,
              help="Conda env activated for GPU child jobs (PMA_CONDA_ENV_GPU), e.g. muagene-gpu.")
@click.option("--gpu-image", default=None,
              help="Machine-local path the GPU .sif is pulled to (default ~/.muagene/images/muagene-gpu.sif).")
@click.option("--gpu-image-uri", default=None,
              help="Pinned registry reference the GPU image is PULLED from, e.g. "
                   "docker://<registry>/muagene-gpu:<tag>. No machine builds the image locally. "
                   "Defaults to machine.config gpu_image_uri.")
@click.option("--singularity-module", default=None,
              help="Module to `module load` before singularity exec (e.g. singularityce/3.11.3). "
                   "Defaults to machine.config singularity_module.")
@click.option("--scratch", default=None,
              help="Optional node-local/fast scratch path to bind into the GPU container "
                   "(exported as PMA_GPU_BIND). The run directory and repo root are always "
                   "bound; this is for extra paths a stage writes outside the run dir.")
@click.option("--env-policy", type=click.Choice(["auto", "manual"]), default="auto", show_default=True,
              help="On a missing/stale env at submit: auto-provision (default) or fail loud with the command.")
def configure_execution(
    config_path: str,
    mode: str,
    confirmed_by_user: bool,
    slurm_partition: str | None,
    slurm_account: str | None,
    resources_scale: float | None,
    conda_env: str | None,
    device: str,
    gpu_partition: str | None,
    gpu_gres: str | None,
    gpu_conda_env: str | None,
    gpu_image: str | None,
    gpu_image_uri: str | None,
    singularity_module: str | None,
    scratch: str | None,
    env_policy: str,
) -> None:
    """Record execution mode and write deliverables/plan/config/site.config + hpc.env."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    paths.ensure()

    params_path = str(paths.parameters_yaml)
    prior_mode = provenance.get_value(params_path, "execution.mode", None)
    prior_confirmed = provenance.get_value(params_path, "execution.user_confirmed", False)

    provenance.set_param(
        params_path,
        "execution.mode", mode,
        source="user", confidence="high",
        rationale=f"Execution backend set via configure-execution --mode {mode}.",
    )

    if mode == "local" and device == "gpu":
        raise click.ClickException(
            "GPU is cluster-only: --device gpu requires --mode slurm (use submit). "
            "Local runs use --mode local with the default --device cpu.")

    # On HPC, gpu routes GPU-capable stages to the GPU partition/gres and container
    # (see workflow/resources.smk _GPU_CAPABLE). Stages that are not GPU-capable
    # always run on CPU regardless of this setting. Local mode is CPU-only.
    provenance.set_param(
        params_path,
        "compute.device", device,
        source="user", confidence="high",
        rationale=f"Compute device set via configure-execution --device {device}.",
    )

    # Explicit, auditable record of whether the USER confirmed this execution mode.
    # `run`/`submit` refuse to launch any compute job until this is true (see
    # `gates.enforce_execution_mode_gate`). Recording the mode alone is not enough —
    # the agent must never silently choose local vs HPC.
    #
    # Re-config semantics: an explicit --confirmed-by-user always confirms. Without
    # it, confirmation is PRESERVED only when the mode is unchanged (e.g. bumping
    # --resources-scale on the same backend) — so resource tweaks don't silently
    # un-confirm a run. A *changed* mode (or one never confirmed) resets to
    # unconfirmed, forcing a fresh user confirmation.
    if confirmed_by_user:
        confirmed = True
        confirmed_rationale = ("User explicitly confirmed the execution mode via "
                               "configure-execution --confirmed-by-user.")
    elif prior_confirmed and prior_mode == mode:
        confirmed = True
        confirmed_rationale = (f"Confirmation preserved across re-config of unchanged "
                               f"mode {mode!r} (no --confirmed-by-user needed for "
                               "resource-only changes).")
    else:
        confirmed = False
        confirmed_rationale = ("Execution mode recorded WITHOUT explicit user "
                               "confirmation (--not-confirmed / default, or mode "
                               "changed); run/submit will refuse compute.")
    provenance.set_param(
        params_path,
        "execution.user_confirmed", confirmed,
        source="user", confidence="high",
        rationale=confirmed_rationale,
    )
    if not confirmed:
        click.echo(
            "NOTE: execution mode recorded but NOT user-confirmed. `run`/`submit` "
            "will refuse to launch any compute job until you confirm local vs HPC "
            "with the user and re-run:\n"
            f"  Processing-MuAgent configure-execution --config {paths.run_yaml} "
            f"--mode {mode} --confirmed-by-user",
            err=True,
        )

    # Machine-level infra knobs default from ~/.muagene/machine.config (written once by
    # Execution-MuAgent `init-machine`), so the operator doesn't re-type manager/module/
    # image/env names per run. Precedence: explicit flag > machine.config > env var.
    mc = hpc.load_machine_config()
    settings: dict[str, str | None] = {
        "slurm_partition": slurm_partition or os.environ.get("PMA_SLURM_PARTITION"),
        "slurm_account": slurm_account or os.environ.get("PMA_SLURM_ACCOUNT"),
        "resources_scale": (
            str(int(resources_scale)) if resources_scale is not None
            else os.environ.get("PMA_RESOURCES_SCALE")
        ),
        "conda_env": (conda_env or mc.get("conda_env") or os.environ.get("PMA_CONDA_ENV")
                      or os.environ.get("CONDA_DEFAULT_ENV")),
        "device": device,
        "slurm_gpu_partition": gpu_partition or os.environ.get("PMA_SLURM_GPU_PARTITION"),
        "slurm_gpu_gres": gpu_gres or os.environ.get("PMA_SLURM_GPU_GRES"),
        "gpu_conda_env": gpu_conda_env or mc.get("gpu_conda_env") or os.environ.get("PMA_CONDA_ENV_GPU"),
        "gpu_image": gpu_image or mc.get("gpu_image") or os.environ.get("PMA_GPU_IMAGE"),
        "gpu_image_uri": gpu_image_uri or mc.get("gpu_image_uri") or os.environ.get("PMA_GPU_IMAGE_URI"),
        "singularity_module": (singularity_module or mc.get("singularity_module")
                               or os.environ.get("PMA_SINGULARITY_MODULE")),
        # Optional extra GPU-container bind (-> PMA_GPU_BIND). Flag > machine.config > env.
        "scratch": scratch or mc.get("scratch") or os.environ.get("PMA_GPU_BIND"),
        # Detected infra: None lets Execution auto-detect; machine.config pins them.
        "env_manager": mc.get("manager"),
        "container_runtime": mc.get("container_runtime"),
        "env_policy": env_policy or mc.get("policy") or "auto",
    }

    if mode == "local":
        click.echo(f"Execution mode: local (device={device}; no hpc.env written).")
        return

    if mode == "slurm" and not settings["slurm_partition"]:
        raise click.ClickException(
            "SLURM mode requires --slurm-partition or PMA_SLURM_PARTITION in the environment.")

    # GPU routing prerequisites (cluster-only; preprocessing stages are CPU-only —
    # _GPU_CAPABLE is empty until the integration subagent adds stages). Fail loud
    # rather than silently submitting with a misconfigured partition/env.
    if device == "gpu":
        click.echo(
            "NOTE: --device gpu prepares cluster GPU routing for the integration "
            "subagent (future). Processing-MuAgent preprocessing is CPU-only.",
            err=True,
        )
        # SLURM GPU is container-only: the job runs inside the PULLED image, so fail
        # loud now if the pinned image reference is missing rather than writing
        # image_uri: null and only discovering it at provision/submit
        # (gpu_image_unavailable).
        if mode == "slurm" and not settings["gpu_image_uri"]:
            raise click.ClickException(
                "SLURM --device gpu requires --gpu-image-uri — a pinned registry reference the GPU "
                "image is PULLED from (e.g. docker://<registry>/muagene-gpu:<tag>) — or gpu_image_uri "
                "in ~/.muagene/machine.config (set once via `Execution-MuAgent init-machine`). "
                "No machine builds the image locally.")
        if mode == "slurm" and not settings["slurm_gpu_gres"]:
            raise click.ClickException(
                "SLURM --device gpu requires --gpu-gres (e.g. gpu:A5000:1). Run `hpc-info` to discover "
                "the GPU partition/gres on this cluster.")

    site_cfg = hpc.write_site_config(paths.site_config, mode=mode, settings=settings)
    out = hpc.write_hpc_env(paths.hpc_env_sh, paths.site_config)
    log_event(run_dir, {"stage": "configure_execution", "event": "configured",
                        "mode": mode, "hpc_env": str(out), "site_config": str(site_cfg)})
    click.echo(f"Execution mode: {mode}")
    click.echo(f"Wrote {site_cfg}")
    click.echo(f"Wrote {out}  (derived from site.config)")
    click.echo("Source this file in your shell before submit/run on the cluster:")
    click.echo(f"  source {out}")


@main.command(name="regenerate-locks")
@click.option("--platform", "platforms", multiple=True, default=("linux-64",), show_default=True,
              help="conda platform(s) to lock for. MuAgene is linux-only; default linux-64.")
def regenerate_locks(platforms: tuple[str, ...]) -> None:
    """Regenerate the CPU conda-lock lockfile from workflow/envs/processing.yaml.

    The YAML is the human source of truth; the committed lock is what actually gets
    installed (solve-free, reproducible). Run this AFTER editing processing.yaml, then
    COMMIT the refreshed lock — `validate-env`/`submit` fail loud (`lock_stale_vs_yaml`)
    when the YAML's content hash no longer matches the lock's recorded `# source-sha256:`.
    Lock generation is a science-authoring act, so it lives in Processing-MuAgent.
    Requires conda-lock: `pip install 'Processing-MuAgent[dev]'`.
    """
    import hashlib
    import shutil
    import subprocess

    import yaml

    man = hpc.load_env_manifest()
    cpu = man.get("cpu") or {}
    yaml_path = hpc.REPO_ROOT / cpu["definition"]
    work = (hpc.REPO_ROOT / cpu["lock"]).parent
    if not yaml_path.exists():
        raise click.ClickException(f"CPU env YAML not found: {yaml_path}")
    # The CPU env is rendered with `conda-lock --kind explicit` — a conda-ONLY format that
    # silently drops any `pip:` subsection. A pip dep here would therefore never reach the
    # lock (so a freshly provisioned env would be missing it and fail validate-env at run
    # time, far from this command). Fail loud instead: every dependency must be a conda
    # package. `- pip` itself (a bare string, for `init-machine`'s editable agent installs)
    # is fine; only a `pip:` mapping is rejected.
    spec = yaml.safe_load(yaml_path.read_text()) or {}
    pip_deps = [d for d in (spec.get("dependencies") or [])
                if isinstance(d, dict) and "pip" in d]
    if pip_deps:
        raise click.ClickException(
            f"{yaml_path.name} has a `pip:` subsection, but the CPU lock is rendered with "
            "`conda-lock --kind explicit` (conda-only) — pip deps would be silently dropped "
            "from the lock and missing from every provisioned env. Move those packages to "
            "conda dependencies (all of MuAgene's are on conda-forge/bioconda).")
    if not shutil.which("conda-lock"):
        raise click.ClickException(
            "conda-lock not found. Install dev deps:  pip install 'Processing-MuAgent[dev]'  "
            "(or: pip install conda-lock).")
    # Stamp the lock with the YAML's content hash; the env preflight compares this (not
    # mtimes — git doesn't preserve those) to detect a lock that drifted from the YAML.
    src_hash = hashlib.sha256(yaml_path.read_bytes()).hexdigest()
    for plat in platforms:
        click.echo(f"conda-lock --kind explicit -p {plat} -f {yaml_path}")
        try:
            subprocess.run(["conda-lock", "--kind", "explicit", "-f", str(yaml_path), "-p", plat],
                           cwd=str(work), check=True)
        except (subprocess.SubprocessError, OSError) as exc:
            raise click.ClickException(f"conda-lock failed for {plat}: {exc}")
        produced = work / f"conda-{plat}.lock"        # conda-lock's default explicit name
        dest = work / f"processing.{plat}.lock"        # the manifest's convention
        if produced.exists():
            produced.replace(dest)
        dest.write_text(f"# source-sha256: {src_hash}\n" + dest.read_text())
        click.echo(f"wrote {dest}")
    click.echo("Lockfile(s) regenerated. Commit them alongside processing.yaml.")


@main.command()
@click.argument("stage")
@click.argument("param_kv")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--rationale", default="User revision")
@click.option("--dry-run", is_flag=True,
              help="Preview the parameter change and exactly which artifacts would be "
                   "deleted (plus the current QC thresholds) — mutate nothing.")
def revise(stage: str, param_kv: str, config_path: str, rationale: str, dry_run: bool) -> None:
    """Update one parameter and reset the stage to awaiting_approval.

    PARAM_KV is key=value, e.g. s1_rna_qc.pct_counts_mt_max=10.0
    """
    stage = _canonical_stage(stage)
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    if "=" not in param_kv:
        raise click.ClickException("param_kv must be key=value")
    key, value = param_kv.split("=", 1)
    # Accept both the short form (min_counts_floor=500) and the full form
    # (s1_rna_qc.min_counts_floor=500). Without this normalisation the key is
    # stored bare in parameters.yaml and effective_params() — which looks for
    # "<stage>.<param>" — never finds it, so the revise has no effect on the
    # plan-review preview or on the real stage at runtime.
    if not key.startswith(f"{stage}."):
        key = f"{stage}.{key}"
    try:
        value_parsed = yaml.safe_load(value)
    except Exception:
        value_parsed = value
    if dry_run:
        _revise_dry_run(run_dir, paths, stage, key, value_parsed)
        return
    provenance.set_param(
        str(paths.parameters_yaml),
        key, value_parsed,
        source="user", confidence="high", rationale=rationale,
    )
    approval.mark_awaiting(run_dir, stage)
    log_event(run_dir, {"stage": stage, "event": "revised", "param": key, "value": value_parsed})
    click.echo(f"Revised {key} = {value_parsed!r}; {_display_stage(stage)} is awaiting_approval.")

    # A revise behaves differently depending on which checkpoint is active.
    if not approval.is_approved(run_dir, "plan_review"):
        # Plan-review checkpoint: the plan is not locked yet. The override only
        # tunes a proposed value, so re-render the plan deliverables (overlay) so
        # what the user reviews equals what will run. No QC stage has executed,
        # so there is nothing downstream to invalidate.
        regenerated = _regenerate_plan_deliverables(run_dir)
        if regenerated:
            log_event(run_dir, {"stage": "plan_review", "event": "plan_deliverables_regenerated",
                                "regenerated": regenerated})
            click.echo(
                "Regenerated plan deliverables so the review reflects this revise "
                f"(overlay): {', '.join(regenerated)}. plan_review stays awaiting_approval."
            )
        return

    # Post-approval (QC-review checkpoint): re-running a QC stage requires
    # clearing its stale downstream artifacts AND the post_qc_review gate
    # outputs, or Snakemake reports "Nothing to be done" and silently skips the
    # re-run. Do this deterministically here so the agent never hand-deletes.
    invalidated = _invalidate_qc_downstream(run_dir, stage)
    if invalidated:
        log_event(run_dir, {"stage": stage, "event": "qc_downstream_invalidated",
                            "deleted": invalidated})
        click.echo(
            f"Invalidated {len(invalidated)} stale downstream/gate artifact(s) so the "
            f"re-run regenerates them (incl. the post_qc_review gate). "
            f"Approve {_display_stage(stage)} (and s3_doublets if S1/S2 changed), then submit."
        )

    # Refresh the param-derived preview layer (S0 *_data_explore figures + qc_explore.json)
    # AND re-render the plan deliverables so both track the revise at the QC-review gate too
    # — the same regeneration the plan_review gate runs. _regenerate_plan_deliverables
    # overlays parameters.yaml on the (still-frozen) preprocessing_plan.json, so the plan
    # files never go stale vs the live overrides; the regenerated qc_review_<run>.md remains
    # the authoritative applied-thresholds record. Cheap: reads the persisted
    # *_qc_metrics.parquet (which survive the post-QC cleanup), no heavy reload. This does
    # NOT re-arm the plan_review gate (it calls the renderers, not the plan-review command).
    regenerated = _regenerate_plan_deliverables(run_dir)
    if regenerated:
        log_event(run_dir, {"stage": stage, "event": "deliverables_regenerated_post_qc",
                            "regenerated": regenerated})
        click.echo(
            "Refreshed the S0 QC-exploration preview + plan deliverables to reflect this "
            f"revise (overlay): {', '.join(regenerated)}. The regenerated qc_review_<run>.md "
            "remains the authoritative applied-thresholds record."
        )


def _stage_states(paths: RunPaths) -> list[tuple[str, str, str]]:
    return _sp.stage_states(paths)


@main.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--watch", is_flag=True,
              help="Poll until a review gate needs approval or manifest completes.")
@click.option("--interval", type=float, default=15.0,
              help="Poll interval in seconds when --watch is set.")
def status(config_path: str, watch: bool, interval: float) -> None:
    """Print per-step pipeline state (S1a–S8 + review gates). With --watch, polls until something changes."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    reporting.run_status(paths, watch=watch, interval=interval)


@main.command(name="hpc-status")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def hpc_status(config_path: str) -> None:
    """Report HPC job health, monitor findings, and per-step pipeline state (one-shot).

    This is Processing-MuAgent's single window onto the Execution-MuAgent supervision
    daemon, which is the sole monitor. It reads only structured JSON
    (latest_snapshot.json + latest_submission.json) — health, silence/tolerance,
    findings, kill_action, and supervisor liveness — and prints once, then exits.

    There is no poll loop here: the daemon does the monitoring. This command drives
    the report-and-repoll rule — after `submit`, report this status, then (while the
    job is still running) re-poll on a non-blocking scheduled wakeup after the seconds
    printed on the `Next check:` line, until monitor.pid is removed or a review gate is
    awaiting approval (`Gate signal present`). Report to the user only when the `State:`
    fingerprint changes.
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    reporting.run_hpc_status(paths)


@main.command(name="supervisor-restart")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--kill-existing/--no-kill-existing", default=True, show_default=True,
              help="Kill any running supervisor before starting the new one.")
def supervisor_restart(config_path: str, kill_existing: bool) -> None:
    """Restart the background supervisor daemon without resubmitting the cluster job.

    Use when the supervisor process died mid-run (crash, OOM, site reboot) but the
    cluster job is still active. Reads latest_submission.json and re-invokes
    resume-monitor as a new daemon (no resubmit).

    The supervisor is the kill-on-hang safety layer. Restarting it restores stall
    detection and auto-cancel protection for the running job.
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    sub_path = paths.run_dir / "internal" / "hpc_monitor" / "latest_submission.json"
    if not sub_path.exists():
        raise click.ClickException(
            "No submission recorded for this run. Use `submit` to start a job."
        )
    if kill_existing:
        killed = hpc.kill_existing_supervisor(run_dir)
        if killed:
            click.echo("Stopped existing supervisor.")
    env = hpc._execution_muagent_env()
    if env is None:
        raise click.ClickException(
            "Execution-MuAgent not found. Install it: pip install -e Execution-MuAgent/"
        )
    cmd = [
        sys.executable, "-m", "execution_muagent.cli", "resume-monitor",
        "--run-dir", str(run_dir),
    ]
    result = hpc.start_supervisor_daemon(run_dir, cmd, env)
    if result is None:
        raise click.ClickException("Failed to start supervisor daemon.")
    pid = result["pid"]
    log = result["log"]
    log_event(run_dir, {
        "stage": "supervisor_restart", "event": "restarted",
        "supervisor_pid": pid, "supervisor_log": log,
    })
    click.echo(f"Supervisor restarted (PID {pid}), logging to {log}")
    click.echo(f"Report status: Processing-MuAgent hpc-status --config {paths.run_yaml}")


@main.command(name="plan-review")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--intro", "intro_text", default=None,
              help="Introductory paragraph to prepend before the Summary section.")
@click.option("--intro-context", "intro_context_only", is_flag=True, default=False,
              help="Print the intro context JSON and exit without writing plan_review.md.")
def plan_review_cmd(config_path: str, intro_text: str | None, intro_context_only: bool) -> None:
    """Render and write the merged plan-review markdown (summary + appendix).

    Also writes per-stage job spec YAMLs to internal/stage_meta/ so Execution-MuAgent
    can read science intent, resource hints, and progress_timeout_hint per stage.

    This command is a *renderer*: it requires the planning compute (P1 → S0)
    to have finished and produced preprocessing_plan.json. Calling it before that
    would emit placeholder deliverables and a false awaiting_approval signal.
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    if intro_context_only:
        import json as _json
        click.echo(_json.dumps(_pr.build_intro_context(run_dir), indent=2))
        return
    # Guard against rendering before planning compute has produced the plan.
    missing = []
    if not paths.preprocessing_plan.exists():
        missing.append(str(paths.preprocessing_plan))
    if not paths.validation_report.exists():
        missing.append(str(paths.validation_report))
    if missing:
        raise click.ClickException(
            "Cannot render plan review — S0 ingest has not finished yet.\n"
            f"Missing: {', '.join(missing)}\n"
            "Wait for the planning job (target plan_review_propose) to complete, "
            "then re-run this command."
        )
    text = _pr.render_merged_markdown(run_dir, intro=intro_text)
    click.echo(text)
    out = _pr.write_summary(run_dir, intro=intro_text)
    click.echo(f"\nWritten: {out}")
    html_out = _pr.write_plan_summary_html(run_dir, intro=intro_text)
    click.echo(f"Written: {html_out}")
    # Write per-stage specs; read workflow_branch from plan if available.
    try:
        import json
        plan_path = RunPaths(run_dir).preprocessing_plan
        branch = "paired"
        if plan_path.exists():
            branch = json.loads(plan_path.read_text()).get("workflow_branch", "paired")
        written = _specs.write_stage_specs(run_dir, branch)
        if written:
            click.echo(f"Wrote {len(written)} stage metadata file(s) to {RunPaths(run_dir).stage_meta_dir}/")
    except Exception:
        pass  # spec writing is best-effort; never block plan-review
    # Arm the plan_review gate only when the plan actually exists. The primary
    # gate-arming path is the plan_review_propose Snakemake rule; this CLI path
    # is a re-render convenience after planning compute has finished.
    approval.mark_awaiting(run_dir, "plan_review")


@main.command(name="marker-gene-check")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option(
    "--force-tsne",
    is_flag=True,
    default=False,
    help="Recompute t-SNE even when a valid cache exists.",
)
@click.option(
    "--plot-only",
    is_flag=True,
    default=False,
    help="Write the figure only; do not refresh QC review reports.",
)
@click.argument("genes", nargs=-1, required=True)
def marker_gene_check_cmd(
    config_path: str,
    force_tsne: bool,
    plot_only: bool,
    genes: tuple[str, ...],
) -> None:
    """Generate before/after marker gene expression plots.

    GENES is one or more gene symbols, e.g. CD3E CD20 EPCAM (matched case-insensitively).

    Uses a cached t-SNE embedding when the cell set is unchanged. By default, QC review
    reports are refreshed automatically after plotting. Pass ``--plot-only`` to skip that.
    """
    from .stages import s1a_ambient as _s1a
    run_dir = _resolve_run_dir(config_path)

    if not genes:
        raise click.UsageError("Provide at least one gene symbol.")

    gene_list = list(genes)
    click.echo(f"Checking marker genes: {', '.join(gene_list)}")
    result = _s1a.run_marker_gene_check(
        run_dir, gene_list, force_tsne=force_tsne, refresh_qc=not plot_only,
    )
    if result["found"]:
        click.echo(f"Plotted: {', '.join(result['found'])}")
    else:
        click.echo("No marker genes found in matrix; figure not written.")
    if result["missing"]:
        click.echo(f"Not found in data: {', '.join(result['missing'])}")
    if result["found"]:
        if plot_only:
            click.echo("Figure written (--plot-only: QC reports unchanged).")
        else:
            click.echo("QC reports refreshed.")


@main.command(name="unlock")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
def unlock_cmd(config_path: str) -> None:
    """Remove stale Snakemake locks for a run after confirming no active process."""
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)
    locks = hpc.snakemake_lock_files(paths.snakemake_workdir)
    if not locks:
        click.echo(f"No Snakemake locks found under {paths.snakemake_workdir}.")
        return
    active = hpc.snakemake_processes_for_workdir(paths.snakemake_workdir)
    if active:
        detail = "\n".join(f"  pid {pid}: {args}" for pid, args in active)
        raise click.ClickException(
            "Refusing to unlock while a local Snakemake process references this workdir:\n"
            f"{detail}"
        )
    # Pass the resolved canonical config (absolute) — snakemake --unlock runs with
    # --directory internal/snakemake, so a relative --config would not resolve there
    # (submit uses paths.run_yaml for the same reason). Keeps `unlock --config <rel>`
    # consistent with every other subcommand.
    _unlock_snakemake(run_dir, paths.run_yaml)
    click.echo(f"Unlocked {paths.snakemake_workdir}")


@main.command(name="run")
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--auto-approve", is_flag=True, help="Auto-approve every checkpoint (noninteractive).")
@click.option("--auto-approve-except", "auto_except", multiple=True,
              help="With --auto-approve, do NOT pre-seed the given stage(s). Repeatable. "
                   "Example: --auto-approve-except qc_review")
@click.option("--no-context", is_flag=True, help="Explicit user choice to proceed without biological context; fields marked status=missing.")
@click.option("--marker-genes", "marker_genes_ack",
              type=click.Choice(["defer", "skip"]), default=None,
              help="With --auto-approve: record an explicit marker-gene decision so "
                   "plan_review can be seeded. 'defer' = check at QC review; "
                   "'skip' = decline. Provide actual genes via `revise` instead to run the check.")
@click.option("--target", default="all")
def run_pipeline(config_path: str, auto_approve: bool, auto_except: tuple[str, ...],
                 no_context: bool, marker_genes_ack: str | None, target: str) -> None:
    """Run the DAG LOCALLY. With --auto-approve, checkpoints are unblocked automatically.

    `run` is local-only: it executes on this machine (local mode) or runs the
    login-node localrules (propose / planning / manifest). All cluster job
    submission and monitoring is owned by Execution-MuAgent via `submit` — there
    is no `run --executor slurm` path.

    Use --auto-approve-except <stage> to keep specific gates honoured (e.g.
    qc_review in headless HPC mode).
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)

    _enforce_execution_mode_gate(run_dir, paths)
    mode = provenance.get_value(str(paths.parameters_yaml), "execution.mode", None)
    if mode == "slurm":
        raise click.ClickException(
            f"execution.mode is {mode!r}, but `run` is local-only. Heavy stages "
            "(starting with S0 ingest) must run on a compute node, never the login "
            "node. Submit instead:\n"
            f"  source {paths.hpc_env_sh}\n"
            f"  Processing-MuAgent submit --config {paths.run_yaml} "
            f"--executor {mode} --target {target}"
        )

    _enforce_context_gate(paths, no_context)

    auto_except = tuple(_canonical_stage(s) for s in auto_except)
    if auto_approve:
        # Pre-seed approval sentinels so snakemake can run the DAG end-to-end in a
        # single invocation; --auto-approve-except keeps the listed stages gated.
        kept = set(auto_except)
        _apply_marker_gene_ack(run_dir, marker_genes_ack)
        _seed_approvals(run_dir, HUMAN_CHECKPOINT_STAGES, note="auto-approved", kept=kept)
        if kept:
            click.echo(f"Auto-approved all stages except: "
                       f"{sorted(_display_stage(s) for s in kept)}. "
                       "Snakemake will stop at those gates.")
    _snakemake(["--configfile", str(paths.run_yaml), target], run_dir)


@main.command()
@click.option("--config", "config_path", required=True, type=click.Path(exists=True))
@click.option("--executor", type=CLUSTER_EXECUTOR_CHOICE, required=True,
              help="Scheduler to submit the head-job to (slurm). "
                   "For local foreground runs use `run` (which is local-only).")
@click.option("--target", default=None,
              help="Override the Snakemake target. Omit to auto-infer the first "
                   "incomplete step (e.g. plan_review_propose for planning, "
                   "post_qc_review_propose, all).")
@click.option("--no-context", is_flag=True,
              help="Explicit user choice to proceed without biological context (planning "
                   "submissions only); fields marked status=missing.")
@click.option("--auto-approve", is_flag=True,
              help="Pre-seed all checkpoint sentinels; head-job runs unattended end-to-end.")
@click.option("--auto-approve-except", "auto_except", multiple=True,
              help="With --auto-approve, keep these gates honoured. Repeatable.")
@click.option("--marker-genes", "marker_genes_ack",
              type=click.Choice(["defer", "skip"]), default=None,
              help="With --auto-approve: record an explicit marker-gene decision so "
                   "plan_review can be seeded. 'defer' = check at QC review; "
                   "'skip' = decline. Provide actual genes via `revise` instead to run the check.")
@click.option("--output", "output_log", type=click.Path(), default=None,
              help="Scheduler output-log path for the head-job (optional).")
@click.option("--unlock-stale-locks", is_flag=True,
              help="If Snakemake locks exist and no local process owns this workdir, "
                   "run snakemake --unlock before submitting.")
@click.option("--watch/--no-watch", default=True,
              help="Start a background supervisor daemon that monitors the cluster job and "
                   "cancels it if it hangs (default: on). The daemon survives SSH disconnect "
                   "unless the site uses KillUserProcesses=yes — use tmux/screen there. "
                   "Returns after job submission is confirmed (≤90 s). "
                   "--no-watch: submit only, NO supervisor daemon started — no stall "
                   "detection, no auto-cancel.")
def submit(config_path: str, executor: str, target: str | None, no_context: bool,
           auto_approve: bool, auto_except: tuple[str, ...],
           marker_genes_ack: str | None,
           output_log: str | None, unlock_stale_locks: bool,
           watch: bool) -> None:
    """Submit the snakemake runner as a SLURM head-job.

    This is the ONLY cluster-execution path: Processing-MuAgent prepares the
    head-job spec + site.config and Execution-MuAgent owns submission and
    monitoring (kill-on-hang, hpc-status). The planning phase targets
    ``plan_review_propose`` (auto-inferred), which pulls P1 → S0 as
    Snakemake dependencies and arms the gate at the end of a single head-job.

    Execution-MuAgent is a hard dependency for cluster submission — it renders the
    submission script, submits the head-job, and owns monitoring. If Execution-MuAgent
    is unavailable, this command fails loudly: there is no manual-submission path.

    The head-job runs on a compute node, activates the project conda env, and
    invokes snakemake with the cluster profile. Snakemake then submits per-stage
    child jobs. The head-job exits when the DAG completes or stops at a missing
    approval gate.

    Typical headless workflow on HPC:

        # Run planning interactively (Phase A), then submit the heavy middle:
        Processing-MuAgent submit --config $CFG --executor slurm \\
                --auto-approve --auto-approve-except post_qc_review

        # After QC review, approve and build the handoff (target auto-inferred):
        Processing-MuAgent approve post_qc_review --config $CFG
        Processing-MuAgent submit --config $CFG --executor slurm

        # Verify the handoff, obtain explicit user confirmation, then submit
        # again; the target now resolves to the unattended S4-S8 finish batch:
        Processing-MuAgent submit --config $CFG --executor slurm
    """
    run_dir = _resolve_run_dir(config_path)
    paths = RunPaths(run_dir)

    # System requirement: confirm execution mode with the user before launching ANY
    # cluster job. Fires before approval seeding, target inference, and the
    # site.config check — so resume submissions (S1+) are gated too, not only S0.
    _enforce_execution_mode_gate(run_dir, paths)

    auto_except = tuple(_canonical_stage(s) for s in auto_except)
    if auto_approve:
        kept = set(auto_except)
        _apply_marker_gene_ack(run_dir, marker_genes_ack)
        _seed_approvals(run_dir, HUMAN_CHECKPOINT_STAGES, note="auto-approved (submit)", kept=kept)
        if kept:
            click.echo(f"Auto-approved all stages except: "
                       f"{sorted(_display_stage(s) for s in kept)}.")
        # Tell the head-job's propose rules not to revoke pre-seeded approvals.
        os.environ["PMA_AUTO_APPROVE"] = "1"

    inferred_target = target is None
    resolved_target = target if target is not None else _infer_submit_target(run_dir)

    # Phase 1 biological-context gate — enforced for planning-phase submissions
    # exactly as `run` does. Resume submissions (S1+) skip it: context was
    # already validated when planning ran. Covers both the canonical auto-inferred
    # target (plan_review_propose) and legacy explicit --target s0_ingest_execute.
    if resolved_target in {"s0_ingest_execute", "plan_review_propose"}:
        _enforce_context_gate(paths, no_context)

    phase_seeded = _prepare_submit_approvals(
        run_dir,
        resolved_target,
        inferred_target=inferred_target,
        auto_approve=auto_approve,
        auto_except=auto_except,
    )

    locks = hpc.snakemake_lock_files(paths.snakemake_workdir)
    if locks:
        active = hpc.snakemake_processes_for_workdir(paths.snakemake_workdir)
        if active:
            detail = "\n".join(f"  pid {pid}: {args}" for pid, args in active)
            raise click.ClickException(
                "Snakemake locks exist and a local Snakemake process still references "
                f"{paths.snakemake_workdir}:\n{detail}"
            )
        lock_list = ", ".join(str(p) for p in locks)
        if not unlock_stale_locks:
            raise click.ClickException(
                "Snakemake lock files already exist for this run, so submitting now "
                "would fail with LockException.\n"
                f"Locks: {lock_list}\n"
                "If no scheduler head/child jobs for this run are active, recover with:\n"
                f"  Processing-MuAgent unlock --config {paths.run_yaml}\n"
                "or resubmit with `--unlock-stale-locks`."
            )
        click.echo(f"Unlocking stale Snakemake locks: {lock_list}")
        _unlock_snakemake(run_dir, paths.run_yaml)

    # Archive the prior run's Snakemake logs so a resubmit's PENDING window is not
    # misread as the previous run's failure by `hpc-status` (stage state is derived
    # from the newest per-rule + main logs). No-op on a first submit.
    archived = hpc.archive_prior_run_logs(paths.snakemake_workdir)
    if archived is not None:
        click.echo(f"Archived previous run logs → {archived}")

    out_path = Path(output_log) if output_log else hpc.head_job_log_path(executor)

    if not paths.site_config.exists():
        raise click.ClickException(
            f"site.config not found at {paths.site_config}. "
            "Run `Processing-MuAgent configure-execution --mode slurm ...` first."
        )

    # Regenerate per-stage specs from the CURRENT code + resources before submitting,
    # so Execution-MuAgent's monitor always verifies the stages that will actually run.
    # Per-stage specs are otherwise only written at plan-review time; a resubmit after
    # a code change (e.g. a renamed stage, or a different PMA_RESOURCES_SCALE) would
    # otherwise run against stale specs. Best-effort — never block submission on this.
    try:
        _specs.write_stage_specs(run_dir, provenance.current_branch(str(paths.parameters_yaml)))
    except Exception:
        pass

    # Write the head-job spec so Execution-MuAgent can render + submit it.
    head_spec_path = _specs.write_head_job_spec(run_dir, resolved_target)

    if watch:
        click.echo("Starting supervision daemon (background)...")

    ea_result = hpc.submit_via_execution_muagent(
        head_spec_path,
        paths.site_config,
        run_dir,
        resolved_target,
        watch=watch,
        kill_on_hang=True,
    )
    if ea_result is None:
        raise click.ClickException(
            "Execution-MuAgent is required for cluster submission but is not available "
            "or returned an error.\n"
            "  Install it:  pip install -e Execution-MuAgent/\n"
            "  Then re-run `Processing-MuAgent submit`."
        )

    import re as _re
    if watch:
        pid = ea_result.get("pid", "?")
        mon_log = ea_result.get("log", "")
        # ea_result["job_id"] is the entry the supervisor appended for THIS
        # submission (or None on timeout). Do NOT fall back to the last manifest
        # entry — that can be a stale job_id from a previous head-job.
        job_id = ea_result.get("job_id")
        if not job_id:
            click.echo(
                "Warning: job ID not yet in execution_manifest.jsonl "
                "(scheduler slow or NFS lag). Check `hpc-status` to confirm submission.",
                err=True,
            )
            job_id = "unknown"
        submitted_log_path = hpc.submitted_log_path(executor, out_path, job_id)
        log_event(run_dir, {
            "stage": "submit", "event": "head_job_submitted",
            "executor": executor, "target": resolved_target,
            "job_id": job_id, "auto_approve": auto_approve,
            "kept_gates": sorted(set(auto_except)),
            "phase_auto_approved": phase_seeded,
            "via_execution_agent": True,
            "supervisor_pid": pid, "supervisor_log": mon_log,
            "head_job_log": str(submitted_log_path),
        })
        click.echo(f"Submitted {executor} head-job: {job_id}")
        click.echo(f"  config:     {paths.run_yaml}")
        click.echo(f"  target:     {resolved_target}")
        if phase_seeded:
            click.echo(f"  phase-auto-approved: {', '.join(phase_seeded)}")
        click.echo(f"  log:        {submitted_log_path}")
        click.echo(f"  supervisor: PID {pid}, log → {mon_log}")
        click.echo(
            "\nThe supervisor daemon is the sole monitor (kill-on-hang) and runs in "
            "the background; it survives SSH disconnect (unless the site uses "
            "KillUserProcesses=yes). Do not run a watch loop — report its status on "
            "demand and act when it signals (job terminal, or a review gate awaiting):\n"
            f"  Report status: Processing-MuAgent hpc-status --config {paths.run_yaml}\n"
            "The daemon signals completion by removing internal/hpc_monitor/monitor.pid; "
            "a review gate shows as 'awaiting_approval'."
        )
    else:
        m = _re.search(r"(?:head-job|job[_-]id)[:\s]+(\S+)", ea_result.get("stdout", ""))
        job_id = m.group(1).strip() if m else ea_result.get("stdout", "").strip().splitlines()[-1]
        submitted_log_path = hpc.submitted_log_path(executor, out_path, job_id)
        log_event(run_dir, {
            "stage": "submit", "event": "head_job_submitted",
            "executor": executor, "target": resolved_target,
            "job_id": job_id, "auto_approve": auto_approve,
            "kept_gates": sorted(set(auto_except)),
            "phase_auto_approved": phase_seeded,
            "via_execution_agent": True,
            "head_job_log": str(submitted_log_path),
        })
        click.echo(f"Submitted {executor} head-job: {job_id} (via Execution-MuAgent)")
        click.echo(f"  config:  {paths.run_yaml}")
        click.echo(f"  target:  {resolved_target}")
        if phase_seeded:
            click.echo(f"  phase-auto-approved: {', '.join(phase_seeded)}")
        click.echo(f"  log:     {submitted_log_path}")
        click.echo(
            "\nNote: --no-watch — no supervisor daemon started. "
            "Stalled/hung jobs will NOT be auto-cancelled.\n"
            f"Attach MuAgene monitoring: Processing-MuAgent supervisor-restart "
            f"--config {paths.run_yaml}"
        )


if __name__ == "__main__":
    main()
