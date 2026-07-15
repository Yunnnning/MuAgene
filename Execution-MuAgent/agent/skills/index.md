# Execution-MuAgent skills

`../system_prompt.md` is always loaded; read the skill below on demand.

## Loading order
1. [`../../AGENT.md`](../../AGENT.md) — identity and scope.
2. [`../system_prompt.md`](../system_prompt.md) — always-loaded policy and hard rules.
3. This index — select the needed procedure.
4. One selected skill — operational procedure.
5. [`../tools.md`](../tools.md) and [`../../../contracts/`](../../../contracts/) — consult
   only for command or machine-contract details.

The public README is operator guidance, not an agent instruction source.

| Skill | Trigger | Commands |
|-------|---------|----------|
| [`workflow.md`](workflow.md) | The run-time lifecycle (validate → render → submit → monitor → report) and the operator-facing env-provisioning commands | `execute-spec`, `resume-monitor`, `init-machine`, `provision-env`, `validate-env`, `doctor` |

**Canonical homes — never restate:** the finding codes you emit →
[`../../../contracts/findings.yaml`](../../../contracts/findings.yaml); every run/machine
state file's writer/reader/lifecycle →
[`../../../contracts/state_model.md`](../../../contracts/state_model.md).
