# MuAgene contracts

Machine-readable single source of truth for everything that crosses an agent
boundary or is handed to a downstream consumer. Both `Processing-MuAgent` and
`Execution-MuAgent` reference these files; prose docs (READMEs, system prompts,
skills) link here and must not restate the shapes or codes.

| File | What it pins | Producer → Consumer |
|------|--------------|---------------------|
| [`findings.yaml`](findings.yaml) | The finding-code registry (`{severity, code, message}`) the monitor emits and Processing acts on, plus env-provisioning codes | Execution → Processing / operator |
| [`state_model.md`](state_model.md) | Every run-state / machine-state file: who writes it, who reads it, its lifecycle | both agents |
| [`post_qc_manifest.schema.json`](post_qc_manifest.schema.json) | The post-QC handoff manifest (`muagene.post_qc_handoff/1`), emitter `executor/stages/qc_handoff.py` | Preprocessing → downstream |

## Conventions
- JSON Schema is **draft 2020-12**. `$id` carries the versioned contract name
  (e.g. `muagene.post_qc_handoff/1`); bump it on a breaking shape change.
- Schemas `require` the keys the emitter always writes and allow additional keys
  (forward-compatible). Enums/`const` mirror the emitter exactly.
- Consistency tests (`tests/test_harness_consistency.py` in each agent) assert the
  live code agrees with these files — e.g. every finding code emitted in
  Execution appears in `findings.yaml`, and a representative manifest validates
  against its schema.

## Status / scope
Authored so far: the handoff manifest schema, the findings registry, and the
state model. Remaining cross-boundary schemas to add under this same directory:
`run_yaml`, `parameters`, `site_config`, `head_job`, `stage_meta`,
`latest_snapshot`, `run_manifest` — each reverse-engineered from its emitter and
covered by a round-trip test.
