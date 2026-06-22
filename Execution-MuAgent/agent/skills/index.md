# Execution-MuAgent skills

`../system_prompt.md` is always loaded; read the skill below on demand.

| Skill | Trigger | Commands |
|-------|---------|----------|
| [`workflow.md`](workflow.md) | The run-time lifecycle (validate → render → submit → monitor → report) and the operator-facing env-provisioning commands | `execute-spec`, `resume-monitor`, `init-machine`, `provision-env`, `validate-env`, `doctor` |

**Canonical homes — never restate:** the finding codes you emit →
[`../../../contracts/findings.yaml`](../../../contracts/findings.yaml); every run/machine
state file's writer/reader/lifecycle →
[`../../../contracts/state_model.md`](../../../contracts/state_model.md).
