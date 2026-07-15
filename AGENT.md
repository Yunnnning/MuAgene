---
name: MuAgene
role: Root orchestrator for the Processing- and Execution-MuAgent instruction sets
registry: muagene.agents.yaml
processing_manifest: Processing-MuAgent/AGENT.md
execution_manifest: Execution-MuAgent/AGENT.md
contracts_dir: contracts
---

# MuAgene

MuAgene coordinates preprocessing workflows for single-cell multi-omics data by orchestrating two specialized agents: Processing-MuAgent and Execution-MuAgent. This repository provides the starting point for agent instructions, contract definitions, and operational boundaries required to deterministically manage all preprocessing stages from data intake through final deliverables. Consult the linked instruction sets for specific procedures according to run state or operational need.

## Agent boundary

- [Processing-MuAgent](Processing-MuAgent/AGENT.md) is user-facing. It owns scientific
  intent, intake, planning, review gates, parameter provenance, and final deliverables.
- [Execution-MuAgent](Execution-MuAgent/AGENT.md) is platform-facing. It validates
  Processing-authored specs, provisions environments, submits and monitors SLURM jobs,
  and reports structured findings. It does not make scientific decisions or contact the
  user during a run.
- The machine-readable registry and handoff direction live in
  [`muagene.agents.yaml`](muagene.agents.yaml). Contract shapes and state ownership live
  in [`contracts/`](contracts/).

Processing is the sole user interface during preprocessing. Cluster delegation is a
CLI subprocess and structured-file handoff, not a second LLM conversation.

## Instruction loading order

1. Read this root file for composition, boundaries, and terminology.
2. Read [`Processing-MuAgent/AGENT.md`](Processing-MuAgent/AGENT.md), then its
   always-loaded [`agent/system_prompt.md`](Processing-MuAgent/agent/system_prompt.md).
3. Enter the Processing
   [`agent/skills/index.md`](Processing-MuAgent/agent/skills/index.md) router and load
   only the skill selected by observable run state.
4. Consult Processing [`agent/tools.md`](Processing-MuAgent/agent/tools.md) only for
   command contracts and [`contracts/`](contracts/) only for shared state or handoff
   details.
5. Load Execution instructions only for platform operation: its
   [`AGENT.md`](Execution-MuAgent/AGENT.md), system prompt, skill index, and selected
   procedure.

Public READMEs are user/operator guidance, not agent policy. If prose conflicts with a
canonical source below, use the canonical source and report the drift.

## Canonical terminology

- **Processing-MuAgent**: the scientific agent role and primary installed CLI name.
- **`executor`**: the supported short CLI alias for `Processing-MuAgent`, and the Python
  package under `Processing-MuAgent/executor/`. In skills, `executor <command>` means a
  shell invocation of that CLI; it does not mean a model-native tool or arbitrary script.
- **Processing skill**: one ordered Markdown procedure under
  `Processing-MuAgent/agent/skills/`. Numeric filename prefixes encode happy-path order;
  frontmatter `name` values remain stable semantic IDs.
- **`executor/stages/*.py`**: scientific stage implementations called by Snakemake rules.
- **Processing `workflow/`**: the Snakemake DAG, rules, environments, and scheduler
  profiles. It is not an agent skill directory.
- **Execution-MuAgent**: the platform CLI; **`execution_muagent`** is its Python package.
- **Execution `agent/skills/workflow.md`**: the platform agent's single operational
  procedure. It is unrelated to Processing's Snakemake `workflow/` directory.
- **tool contract**: a Markdown description in `agent/tools.md` of a live Click command.
  MuAgene has no in-repository JSON function-call or MCP tool registry.

## Runtime path

The agent reads one Processing skill, which names allowed CLI commands in
`calls_tools`. The external host runs the selected command through a shell:

```text
Processing skill
  -> executor CLI
  -> local: Snakemake -> executor/stages/*.py
  -> SLURM: Execution-MuAgent execute-spec -> head job -> Snakemake -> child jobs
  -> Execution latest_snapshot.json -> executor hpc-status -> Processing -> user
```

Skills describe when to act; CLIs and Snakemake perform and constrain the action. Never
hand-edit run state. Use `executor` or `Execution-MuAgent` commands according to
[`contracts/state_model.md`](contracts/state_model.md).

## Source-of-truth precedence

- Agent composition and vocabulary: this file.
- Agent identity and scope: each component `AGENT.md`.
- Behavioral hard rules: each component `agent/system_prompt.md`.
- Current procedure: the selected skill reached through `agent/skills/index.md`.
- CLI flags and behavior: live Click implementations; `agent/tools.md` is their contract.
- Stage topology and branch membership: `Processing-MuAgent/executor/pipeline.py`.
- QC defaults: `Processing-MuAgent/executor/defaults.py`.
- Paths: `Processing-MuAgent/executor/run_paths.py`.
- Finding codes and state ownership: `contracts/findings.yaml` and
  `contracts/state_model.md`.
- Cross-agent registry and handoff versions: `muagene.agents.yaml` and schemas under
  `contracts/`.

Do not duplicate canonical values or procedures in a higher-level file. Link to their
owner instead.
