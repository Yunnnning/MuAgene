"""parameters.yaml reader/writer with provenance rules.

Schema: flat map keyed by <stage>.<param>. Each entry has:
  value, source (user|derived|default|recommended|inferred),
  confidence (high|medium|low), rationale, assumptions, method? (required for derived/inferred),
  approved_by?, approved_at?, revision_of?
"""
from __future__ import annotations

import datetime as _dt
from pathlib import Path
from typing import Any

import yaml

ALLOWED_SOURCES = {"user", "derived", "default", "recommended", "inferred"}
ALLOWED_CONFIDENCE = {"high", "medium", "low"}


def _utcnow() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load(path: Path | str) -> dict[str, Any]:
    p = Path(path)
    if not p.exists():
        return {}
    with p.open() as f:
        return yaml.safe_load(f) or {}


def save(path: Path | str, params: dict[str, Any]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as f:
        yaml.safe_dump(params, f, sort_keys=True, allow_unicode=True)


def validate_entry(key: str, entry: dict[str, Any]) -> None:
    required = {"value", "source", "confidence", "rationale"}
    missing = required - entry.keys()
    if missing:
        raise ValueError(f"param {key}: missing required fields {missing}")
    if entry["source"] not in ALLOWED_SOURCES:
        raise ValueError(f"param {key}: source={entry['source']!r} not in {ALLOWED_SOURCES}")
    if entry["confidence"] not in ALLOWED_CONFIDENCE:
        raise ValueError(f"param {key}: confidence={entry['confidence']!r} not in {ALLOWED_CONFIDENCE}")
    if entry["source"] in {"derived", "inferred"}:
        method = entry.get("method") or {}
        if not method.get("name"):
            raise ValueError(f"param {key}: source={entry['source']} requires method.name")
    if entry["source"] == "user" and "method" in entry:
        raise ValueError(f"param {key}: source=user forbids method")


def set_param(
    path: Path | str,
    key: str,
    value: Any,
    *,
    source: str,
    confidence: str,
    rationale: str,
    assumptions: list[str] | None = None,
    method: dict[str, Any] | None = None,
    approved_by: str | None = None,
) -> None:
    params = load(path)
    prior = params.get(key)
    entry: dict[str, Any] = {
        "value": value,
        "source": source,
        "confidence": confidence,
        "rationale": rationale,
        "assumptions": assumptions or [],
    }
    if method is not None:
        entry["method"] = method
    if approved_by is not None:
        entry["approved_by"] = approved_by
        entry["approved_at"] = _utcnow()
    if prior is not None and prior.get("value") != value:
        entry["revision_of"] = {k: prior.get(k) for k in ("value", "source", "approved_at")}
    validate_entry(key, entry)
    params[key] = entry
    save(path, params)


def current_branch(params_path: Path | str, default: str = "paired") -> str:
    """Return the effective workflow_branch for a run.

    Order of precedence:
      1. `plan.workflow_branch` — committed by S0 after detection (+ optional
         confirmation of a prior declaration).
      2. `plan.workflow_branch_declared` — user assertion written by
         `executor declare-branch` before S0 runs. Used by Snakemake input
         functions during dry-run (before S0 has committed the final value).
      3. `default` — typically "paired" to preserve legacy behaviour.
    """
    committed = get_value(params_path, "plan.workflow_branch", None)
    if committed:
        return committed
    declared = get_value(params_path, "plan.workflow_branch_declared", None)
    if declared:
        return declared
    return default


def get_value(path: Path | str, key: str, default: Any = None) -> Any:
    params = load(path)
    if key in params:
        return params[key]["value"]
    return default


def effective_value(
    path: Path | str,
    plan_params: dict[str, Any],
    stage: str,
    name: str,
    default: Any = None,
) -> Any:
    """Effective value for ``<stage>.<name>``.

    The single overlay rule shared by the QC stages and the plan-review
    renderer: a user ``revise`` recorded in parameters.yaml wins over the frozen
    plan default. Falls back to the plan entry's value, then ``default``.
    """
    v = get_value(path, f"{stage}.{name}", None)
    if v is not None:
        return v
    entry = plan_params.get(name)
    if isinstance(entry, dict) and "value" in entry:
        return entry["value"]
    return default


def effective_params(
    path: Path | str, plan_params: dict[str, Any], stage: str
) -> dict[str, Any]:
    """``plan_params`` overlaid with parameters.yaml overrides → ``{name: entry}``.

    Each parameter's entry is replaced by the parameters.yaml entry when one
    exists for ``<stage>.<name>`` (so a ``revise`` is reflected). Keys present
    only in parameters.yaml are included too — harmless, since consumers look up
    by the plan's parameter names. Use when a whole stage's parameter set is
    consumed (e.g. the QC-exploration preview).
    """
    eff = dict(plan_params)
    prefix = f"{stage}."
    for key, entry in load(path).items():
        if key.startswith(prefix) and isinstance(entry, dict) and "value" in entry:
            eff[key[len(prefix):]] = entry
    return eff
