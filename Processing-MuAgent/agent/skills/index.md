# Processing-MuAgent skills — router

`../system_prompt.md` is always loaded and points here. This file is the **router**:
identify the current stage from observable state (`executor status` / which gate is
`awaiting_approval` / whether a cluster job is running), then read **only that stage's
skill**. Don't load the whole flow at once — one stage skill per turn.
Per-command tool contracts: [`../tools.md`](../tools.md).

## State → skill router

Pick the first row whose condition matches the current observable state.

| Observable state | Read skill | Domain |
|---|---|---|
| New conversation / no run dir yet | [`entry_declare.md`](entry_declare.md) | entry |
| Analysis type + run_dir known; `run.yaml` not written yet | [`inputs_intake.md`](inputs_intake.md) | intake |
| Run scaffolded (`init` done), branch declared, exec-mode confirmed; `plan_review` not approved | [`plan_confirm.md`](plan_confirm.md) | plan |
| `plan_review` approved; QC stages (S1a–S3) not yet complete | [`run_execution.md`](run_execution.md) | execute |
| `status`: `post_qc_review` is `awaiting_approval` | [`qc_review_and_revise.md`](qc_review_and_revise.md) | QC |
| `post_qc_review` approved; finish batch (S4–S8) running | [`downstream_dimred_clustering.md`](downstream_dimred_clustering.md) | downstream |
| `manifest` complete | [`completion_handoff.md`](completion_handoff.md) | finish |
| Any cluster job running (during any batch) | [`hpc_monitoring.md`](hpc_monitoring.md) | monitor |
| Any finding/error in the snapshot, or a raised executor error | [`troubleshooting.md`](troubleshooting.md) | recovery |

The happy path is linear (top to bottom). `hpc_monitoring` and `troubleshooting` are
**cross-cutting** — enter them from any compute stage and return to the row matching the
new state. **After every gate approval, re-run `executor status` and re-enter this router.**

**QC threshold revision is gate-scoped, not its own row.** Before `plan_review` is approved,
a "change the thresholds" request is handled inside [`plan_confirm.md`](plan_confirm.md)
(non-destructive there — it just re-renders the plan). [`qc_review_and_revise.md`](qc_review_and_revise.md)
is entered **only** when `post_qc_review` is `awaiting_approval`; it owns the post-run gate and
is the canonical reference for the revise keys both gates link to.

## Skill frontmatter contract

Every skill opens with this YAML block — its machine-readable contract. Read the
frontmatter to confirm you are in the right stage; read the body only when you act.

```yaml
---
name: <skill_id>
domain: <entry|intake|plan|execute|QC|monitor|downstream|finish|recovery>
purpose: <one line>
activation: <observable state that selects this skill>      # when to enter
inputs:  [<state files / prior outputs consumed>]           # input contract
outputs: [<state files this skill produces>]                # output contract
calls_tools:     [<executor subcommands>]
reads_contracts: [<names under ../../contracts/>]
writes_state:    [<state files mutated — via the CLI only>]
handoff: { next: <skill|STOP>, when: <advance condition>, on_error: troubleshooting }
---
```

## Canonical homes — never restate these elsewhere; link instead

| Fact | Single home |
|---|---|
| QC default values / fixed Leiden resolutions | `executor/defaults.py` |
| Finding codes, state-file lifecycle, handoff schemas | [`../../../contracts/`](../../../contracts/) |
| Marker-gene "never invent genes" rule | [`qc_review_and_revise.md`](qc_review_and_revise.md) |
| QC revise reference (keys, `*_override`, skip recipes, binding-constraint diagnosis) — used at both gates | [`qc_review_and_revise.md`](qc_review_and_revise.md) |
| Report-and-repoll monitoring rule | [`hpc_monitoring.md`](hpc_monitoring.md) |
| Execution-mode intake heuristics (file-size → scale) | [`inputs_intake.md`](inputs_intake.md) |
| Error → remedy scenarios | [`troubleshooting.md`](troubleshooting.md) |
