# Execution-MuAgent

Execution-MuAgent is MuAgene's platform execution layer. It prepares and supervises
SLURM jobs and provisions the environments they need. It is deliberately science-free:
[Processing-MuAgent](../Processing-MuAgent/README.md) decides what preprocessing to run
and remains the user-facing agent.

## Responsibilities

Execution-MuAgent owns everything from a confirmed cluster specification to a monitored
job, plus machine-level environment provisioning:

- bootstraps and validates the MuAgene runtime environment;
- validates Processing-authored job specifications;
- renders and submits SLURM jobs;
- monitors scheduler and workflow progress;
- verifies declared outputs;
- reports structured findings to Processing-MuAgent;
- cancels work only after evidence confirms an unhealthy or filesystem-hung job.

It does not choose scientific methods, change preprocessing parameters, modify job
specifications, talk to users during a run, or resubmit failed work.

### Who runs which commands?

| Goal                                   | User-facing Interface      | Actual Execution (Owner)        |
|-----------------------------------------|---------------------------|---------------------------------|
| Set up a new machine                   | Operator (direct CLI)     | Execution-MuAgent               |
| Provision or validate environments      | Operator (direct CLI)     | Execution-MuAgent               |
| Submit a preprocessing run              | Processing-MuAgent        | Execution-MuAgent (behind scenes)|
| Check run status                        | Processing-MuAgent        | Execution-MuAgent (reports back)|

> **Note:**  The user interacts exclusively with Processing-MuAgent during regular preprocessing runs. However, Execution-MuAgent is responsible for actually rendering, submitting, and supervising jobs on the cluster as instructed by Processing-MuAgent. Direct operator interaction with Execution-MuAgent is limited to machine setup and environment management.  

## Requirements

- Linux
- SLURM access for cluster execution
- `micromamba`, `mamba`, or `conda`
- Python 3.10–3.12
- Processing-MuAgent and Execution-MuAgent checked out as sibling directories


## Installation

From the Execution-MuAgent repository:

```bash
bash scripts/bootstrap.sh --processing-repo ../Processing-MuAgent
conda activate muagene
```

The bootstrap script:

1. detects the available environment manager and relevant platform capabilities.
2. creates or updates the integrated `muagene` environment.
3. installs both MuAgene agent packages.
4. validates the installation.
5. records reusable machine settings.

CPU setup is the default. CPU environments are created directly from the committed environment lock file; no dependency-solving or package downloads are required on the target machine. 

GPU environments use a pre-built container image pulled from a registry (not built locally). To prepare GPU infrastructure as well, supply a pinned GPU image reference and any required container module:

```bash
bash scripts/bootstrap.sh \
  --processing-repo ../Processing-MuAgent \
  --device both \
  --gpu-image-uri docker://REGISTRY/IMAGE:TAG \
  --singularity-module MODULE
```

> **Note:** All current preprocessing steps run exclusively on CPU. GPU workflows are not yet supported; GPU configuration exists only in preparation for future features (such as multiomic integration) and is not used in present workflows.

## Operator commands

| Command | Purpose |
|---|---|
| `init-machine` | Probe, provision, install, and validate a machine |
| `provision-env` | Create or refresh CPU and/or GPU environments |
| `validate-env` | Check environment identity and required imports |
| `doctor` | Report platform capabilities and environment health |

These commands are operator-facing and fail loudly with actionable diagnostics. The
`report` command is an advanced read-only debug helper; normal run status comes through
Processing-MuAgent.

## What happens during a cluster run?

When Processing-MuAgent submits work, Execution-MuAgent:

1. validates the requested resources, inputs, platform, and environments;
2. renders a scheduler script and submits it to SLURM;
3. records the accepted job and starts background supervision;
4. watches scheduler state, workflow progress, and declared outputs;
5. returns structured progress or failure findings to Processing-MuAgent.

Policy rejections are not retried. Transient submission failures receive a limited retry.
A quiet job is investigated before any cancellation. Execution-MuAgent never resubmits;
Processing-MuAgent presents the evidence and recovery choice to the user.

## Troubleshooting

- **Machine setup fails:** run `Execution-MuAgent doctor`, correct the reported platform
  or environment issue, then rerun bootstrap or provisioning.
- **An environment is missing or stale:** run `Execution-MuAgent provision-env` for the
  affected device, or let the configured automatic policy reconcile it on submission.
- **The supervisor stopped while the SLURM job is still active:** use
  Processing-MuAgent's supervisor-restart action. It resumes monitoring without
  submitting a second job.
- **Your cluster kills background processes at logout:** run the submission session
  inside `tmux` or `screen`.
- **A run fails:** inspect status through Processing-MuAgent. Execution-MuAgent reports
  evidence but does not choose or apply the scientific recovery.

## Project context

The public README intentionally omits monitor algorithms and internal state layouts.
Canonical details live in:

- [AGENT.md](AGENT.md) — concise role and contract;
- [agent instructions](agent/skills/index.md) — runtime procedures;
- [tool contracts](agent/tools.md) — command behavior;
- [shared contracts](../contracts/) — state ownership, schemas, and finding codes;
- [Processing-MuAgent](../Processing-MuAgent/README.md) — the user-facing preprocessing
  workflow.

Execution-MuAgent is the infrastructure component of [MuAgene](../README.md).
