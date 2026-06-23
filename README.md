# MuAgene

A two-subagent framework for **single-cell multi-omics preprocessing**. One agent owns the
science; the other owns the machine. They communicate through versioned, machine-readable
contracts — not prose.

## The two agents

| Agent | Owns | Manifest |
|-------|------|----------|
| **[Processing-MuAgent](Processing-MuAgent/AGENT.md)** | Scientific intent — branch selection, biological context, the QC/preprocessing plan, parameter provenance, the two human-in-the-loop gates, deterministic deliverables. Pipeline `P1 → S0 → [plan_review] → S1a..S3 → [post_qc_review] → s_handoff + S4..S8 → manifest`. | `Processing-MuAgent/AGENT.md` |
| **[Execution-MuAgent](Execution-MuAgent/AGENT.md)** | Platform mechanics + machine infra — validate spec → render → submit (SLURM) → monitor → report; plus environment provisioning (CPU conda-lock, pull-only GPU container). Science-free; never contacts the user during a run. | `Execution-MuAgent/AGENT.md` |

The boundary is declared in [`muagene.agents.yaml`](muagene.agents.yaml): Processing writes
`site.config` + the per-stage specs; Execution writes `latest_snapshot.json` (findings +
monitor state + kill action) that Processing reads and acts on.

## Repository map

```
MuAgene/
├── README.md                 # you are here
├── muagene.agents.yaml       # agent registry + the inter-agent contract boundary
├── contracts/                # SHARED single source of truth (machine-readable)
│   ├── findings.yaml         #   inter-agent finding-code registry
│   ├── state_model.md        #   every run/machine state file: writer / reader / lifecycle
│   └── *.schema.json         #   cross-boundary + handoff artifact schemas
├── Processing-MuAgent/
│   ├── AGENT.md              # manifest (identity, scope, I/O, hard rules)
│   ├── agent/                # system_prompt.md + skills/ (on-demand procedures) + tools.md
│   ├── executor/             # the Click CLI + science stages (executor/defaults.py = QC-default SSOT)
│   └── workflow/             # Snakemake DAG + envs
└── Execution-MuAgent/
    ├── AGENT.md
    ├── agent/                # system_prompt.md + skills/ + tools.md
    └── execution_muagent/    # the orchestration/monitoring CLI
```

## How the harness is organized
Each concern has exactly one home, referenced everywhere else (no restatement):
- **Identity & policy** → each agent's `AGENT.md` + slim `agent/system_prompt.md`.
- **Procedures** → `agent/skills/` (progressive disclosure; start at `skills/index.md`).
- **Tool behavior** → `agent/tools.md` (per CLI command).
- **Cross-boundary shapes, finding codes, state lifecycle** → `contracts/`.
- **QC default values** → `Processing-MuAgent/executor/defaults.py`.

Drift is caught by `tests/test_harness_consistency.py` in each agent (e.g. plan defaults ==
`defaults.py`, every emitted finding code is registered, the handoff manifest validates
against its schema).

## Quickstart

**Bootstrap a machine once** (operator-facing; creates the `muagene` env + installs both packages):
```bash
Execution-MuAgent init-machine --processing-repo /path/to/Processing-MuAgent --device cpu
```

**Run a preprocessing job** (the agent drives this conversationally; the CLI underneath):
```bash
executor init --config run.yaml                 # scaffold the run dir
executor declare-branch paired --config $CFG    # paired | separate | rna_only | atac_only
executor configure-execution --config $CFG --mode local --confirmed-by-user   # or --mode slurm
executor plan-review --config $CFG              # gate #1 — review deliverables/plan/plan_review_<run>.md
executor approve plan_review --config $CFG
executor run --config $CFG                      # local; or: executor submit (cluster, via Execution-MuAgent)
#   ... at gate #2 (post_qc_review): review, optionally `revise`, then approve ...
```

See each agent's `agent/system_prompt.md` and `agent/skills/index.md` for the full flow.

## Conventions
- All run state mutates **only** through the `executor` / `Execution-MuAgent` CLIs.
- Raw inputs stay pristine; derived files are written alongside, never overwritten.
- Failures fail **loud** (no silent degradation); environment/execution errors are fixed at
  the environment level and documented.
- Tests: per agent, `PYTHONPATH=. conda run -n muagene python -m pytest tests/ -q`.
