# Processing-MuAgent skills

Progressive-disclosure procedures. `../system_prompt.md` is always loaded; read the skill
for the current phase **on demand**. Each skill names its trigger and the `executor`
commands it calls (per-command contracts in [`../tools.md`](../tools.md)).

| Skill | Trigger | Calls (executor) |
|-------|---------|------------------|
| [`entry_declare.md`](entry_declare.md) | New conversation — pick the analysis type | — (dialogue only) |
| [`inputs_intake.md`](inputs_intake.md) | Type known — collect paths + biological context + execution mode (local/HPC, device) | `init`, `declare-branch`, `hpc-info`, `configure-execution` |
| [`workflow.md`](workflow.md) | The end-to-end flow: confirm the plan, then run with checkpoints | `plan-review`, `approve`, `revise`, `run`/`submit`, `status` |
| [`qc_review_and_revise.md`](qc_review_and_revise.md) | `post_qc_review` gate, or any QC-threshold change | `status`, `revise`, `approve`, `marker-gene-check` |
| [`hpc_monitoring.md`](hpc_monitoring.md) | After `submit` — track cluster job health | `hpc-status` (report-and-repoll) |

**Canonical homes — never restate these elsewhere:** QC default values →
`executor/defaults.py`; finding codes, state-file lifecycle, and handoff schemas →
[`../../../contracts/`](../../../contracts/); the marker-gene "never invent genes" rule →
[`qc_review_and_revise.md`](qc_review_and_revise.md); the report-and-repoll monitoring rule →
[`hpc_monitoring.md`](hpc_monitoring.md).
