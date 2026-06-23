---
name: hpc_monitoring
domain: monitor
purpose: Canonical report-and-repoll rule — report one-shot hpc-status, re-poll on a non-blocking scheduled wakeup, drive the gate on signal. Other skills link here.
activation: after executor submit, while a cluster job is running
inputs: [internal/hpc_monitor/latest_snapshot.json, internal/hpc_monitor/monitor.pid]
outputs: [user-facing status report]
calls_tools: [hpc-status]
reads_contracts: [latest_snapshot, findings]
writes_state: []
handoff: { next: re-enter index.md router, when: monitor.pid gone or gate awaiting_approval, on_error: troubleshooting }
---

# Skill: hpc_monitoring — report-and-repoll

**Trigger:** after `executor submit`, while a cluster job is running. This is the
canonical statement of the monitoring rule; other skills link here instead of restating it.

Job monitoring is owned by the **Execution-MuAgent supervision daemon** that `submit`
starts in the background (it refreshes `internal/hpc_monitor/latest_snapshot.json` each
poll and owns kill-on-hang). Do **not** run a blocking `executor hpc-status --watch` or
`tail -f` — that duplicates the daemon and blocks the session.

## Procedure
1. Run one-shot `executor hpc-status --config $CFG`; relay the status to the user.
2. Read the snapshot's `Next check:` line and arm a **non-blocking** scheduled wakeup at
   that cadence (the daemon heartbeat + a small buffer). Never a foreground loop.
3. On wake, re-poll; re-report **only** when the `State:` fingerprint line changed.
4. Stop polling and drive the gate when `internal/hpc_monitor/monitor.pid` is gone **or**
   `hpc-status` prints `Gate signal present`.

## The contract
`latest_snapshot.json` is the machine contract Processing consumes — `findings`,
`monitor_state`, `kill_action`, `error_context`, and the `interval_s` /
`next_recheck_after_s` cadence. Finding codes and their meaning + your recovery action:
[`../../../contracts/findings.yaml`](../../../contracts/findings.yaml). Full state-file
lifecycle: [`../../../contracts/state_model.md`](../../../contracts/state_model.md).
Processing owns all recovery (escalate / fix / resubmit); Execution never contacts the user
during a run. `latest_report.md` is daemon-internal — never parse it or show it to the user.
