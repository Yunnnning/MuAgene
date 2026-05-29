"""Build the S7 resolution-review deliverable: a notebook + a static HTML report.

Invoked by the `s7_clustering_propose` rule after the sweep + adjacency report are
written. Produces:
    <run_dir>/deliverables/checkpoint/resolution_review/resolution_review.ipynb
    <run_dir>/deliverables/checkpoint/resolution_review/resolution_review.html

The HTML is rendered statically from the resolution_summary.md + sweep table +
adjacency report — no notebook execution required, so it's safe to call from the
local propose rule on a login node. The .ipynb is provided for power users who
want to re-run with custom resolutions interactively (via JupyterLab or papermill).

Design notes:
- The HTML must be openable in any browser with no tooling; it is the primary
  artifact for HPC users without a Jupyter session.
- The notebook is self-contained: RUN_DIR is baked in at generation time and
  can be overridden at runtime via the PMA_RUN_DIR env var.
- Avoid executing the notebook here. Re-running Leiden at adjacent resolutions
  is the user's choice (interactively in Jupyter, or via
  `Processing-MuAgent resolution-compare`); we just provide the scaffold.
"""
from __future__ import annotations

import hashlib
import html
import json
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Notebook cells (nbformat v4.5 — no nbformat dependency)
# ---------------------------------------------------------------------------

def _cell_id(source: str, salt: str) -> str:
    return hashlib.sha1((salt + "::" + source).encode()).hexdigest()[:12]


def _md(source: str, salt: str = "") -> dict[str, Any]:
    return {"cell_type": "markdown",
            "id": _cell_id(source, salt or f"md:{source[:20]}"),
            "metadata": {},
            "source": source.splitlines(keepends=True)}


def _code(source: str, salt: str = "") -> dict[str, Any]:
    return {"cell_type": "code",
            "id": _cell_id(source, salt or f"code:{source[:20]}"),
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": source.splitlines(keepends=True)}


_CELL_SETUP = """\
import os
import json
from pathlib import Path
import pandas as pd
import yaml
from IPython.display import Markdown, display

RUN_DIR = Path(os.environ.get("PMA_RUN_DIR", "__BAKED_RUN_DIR__"))
ART = RUN_DIR / "internal" / "artifacts" / "s7_clustering"
SUMMARY_MD = RUN_DIR / "deliverables" / "checkpoint" / "resolution_review" / "resolution_summary.md"
PARAMS_YAML = RUN_DIR / "internal" / "parameters.yaml"

sweep = pd.read_parquet(ART / "sweep.parquet")
params = yaml.safe_load(PARAMS_YAML.read_text()) or {}
adjacency = json.loads((ART / "adjacency_report.json").read_text())

def _get(key):
    rec = params.get(key)
    return rec.get("value") if isinstance(rec, dict) else None

rna_res = _get("s7_clustering.rna.resolution")
atac_res = _get("s7_clustering.atac.resolution")

print(f"Recommended resolutions — RNA: {rna_res}, ATAC: {atac_res}")
print(f"Run dir: {RUN_DIR}")
"""


_CELL_SHOW_SUMMARY = """\
display(Markdown(SUMMARY_MD.read_text()))
"""


_CELL_SWEEP_TABLE = """\
for modality in ("rna", "atac"):
    sub = sweep[sweep["modality"] == modality]
    if sub.empty:
        continue
    chosen = rna_res if modality == "rna" else atac_res
    print(f"=== {modality.upper()} sweep (chosen: {chosen}) ===")
    display(sub[["resolution", "n_clusters", "silhouette", "seed_stability_ari"]]
            .reset_index(drop=True))
"""


_CELL_COMPARE_SIDE_BY_SIDE = """\
# Re-cluster at the recommended ± one user-chosen value and render side-by-side
# UMAPs. Adjust the tuples below to explore different alternatives.
#
# This re-runs Leiden but does NOT modify the approved cluster labels on disk.
from executor import resolution_compare as _rcmp

# Default to comparing recommended vs +0.2 step; adjust to taste.
RNA_PAIR  = (float(rna_res or 1.0),  float(rna_res or 1.0) + 0.2)
ATAC_PAIR = (float(atac_res or 0.8), float(atac_res or 0.8) + 0.2)

rna_report  = _rcmp.compare_rna(RUN_DIR, RNA_PAIR)
atac_report = _rcmp.compare_atac(RUN_DIR, ATAC_PAIR)

for fig in rna_report.get("figures", []) + atac_report.get("figures", []):
    from IPython.display import Image, display
    display(Image(filename=str(fig)))
"""


_CELL_APPROVAL = """\
# To approve the recommended resolutions:
#   Processing-MuAgent approve s7_clustering --config <run_dir>/deliverables/pre_run/config/run.yaml
#
# To revise one (example: bump RNA to 1.2):
#   Processing-MuAgent revise s7_clustering s7_clustering.rna.resolution=1.2 \\
#                              --config <run_dir>/deliverables/pre_run/config/run.yaml
#
# After approval, resubmit:
#   Processing-MuAgent submit --executor pbs   # or local / slurm
print("See the markdown cell above for the exact CLI commands.")
"""


def _resolution_checkpoint_note(run_dir: Path | str) -> str:
    from .. import provenance as _prov
    from ..run_paths import RunPaths
    run_dir = Path(run_dir)
    branch = _prov.get_value(str(RunPaths(run_dir).parameters_yaml),
                             "plan.workflow_branch", "paired") or "paired"
    if branch == "paired":
        return (
            "**Paired multiome:** resolutions are **diagnostic** per-modality Leiden "
            "labels on the joint cell set — they colour UMAPs in `processed.h5mu` but "
            "are not joint integrated clustering."
        )
    if branch == "separate":
        return (
            "**Separate branch:** resolutions set **final** `leiden_rna` / `leiden_atac` "
            "labels in the processed h5ad outputs."
        )
    return (
        f"**{branch} branch:** resolution sets **final** cluster labels in the "
        "processed output."
    )


def _build_cells(run_dir: str) -> list[dict[str, Any]]:
    def bake(s: str) -> str:
        return s.replace("__BAKED_RUN_DIR__", run_dir)

    branch_note = _resolution_checkpoint_note(run_dir)
    return [
        _md(f"# Clustering resolution review\n\n"
            f"Run directory: `{run_dir}`\n\n"
            f"{branch_note}\n\n"
            f"Use the sweep tables below to see how resolution affects n_clusters, "
            f"silhouette, and stability ARI. Approve or revise before S8. Static HTML "
            f"companion: `resolution_review.html`."),
        _md("## Setup"),
        _code(bake(_CELL_SETUP)),
        _md("## Recommendation summary\n\nRendered from the canonical "
            "`deliverables/checkpoint/resolution_review/resolution_summary.md`."),
        _code(_CELL_SHOW_SUMMARY),
        _md("## Sweep tables\n\nFull per-resolution results from `sweep.parquet`."),
        _code(_CELL_SWEEP_TABLE),
        _md("## Side-by-side UMAP comparison (optional)\n\n"
            "Re-clusters at the recommended resolution vs one alternative. "
            "Edit `RNA_PAIR` / `ATAC_PAIR` to compare different values; this "
            "does not modify the approved labels."),
        _code(_CELL_COMPARE_SIDE_BY_SIDE),
        _md("## Approval"),
        _code(_CELL_APPROVAL),
    ]


def build_notebook(run_dir: Path | str) -> dict[str, Any]:
    return {
        "cells": _build_cells(str(run_dir)),
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "mimetype": "text/x-python",
                              "file_extension": ".py", "pygments_lexer": "ipython3",
                              "nbconvert_exporter": "python", "version": "3.10",
                              "codemirror_mode": {"name": "ipython", "version": 3}},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


# ---------------------------------------------------------------------------
# Static HTML rendering (no notebook execution required)
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = """\
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Clustering resolution review</title>
  <style>
    body {{ font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
            max-width: 960px; margin: 2em auto; padding: 0 1em;
            color: #222; line-height: 1.5; }}
    h1, h2, h3 {{ color: #1a3a5c; }}
    code, pre {{ font-family: ui-monospace, "SF Mono", Consolas, monospace;
                 background: #f4f6f8; padding: 0.1em 0.3em; border-radius: 3px; }}
    pre {{ padding: 0.8em; overflow-x: auto; }}
    table {{ border-collapse: collapse; margin: 1em 0; }}
    th, td {{ border: 1px solid #cdd6df; padding: 0.4em 0.7em; }}
    th {{ background: #eef2f7; text-align: left; }}
    .recommended {{ background: #fffceb; font-weight: 600; }}
    .meta {{ color: #666; font-size: 0.92em; margin-bottom: 1.5em; }}
    .alert {{ background: #fff4e5; border-left: 4px solid #f0ad4e;
              padding: 0.6em 0.9em; margin: 1em 0; border-radius: 3px; }}
  </style>
</head>
<body>
<h1>Clustering resolution review</h1>
<p class="meta">Run: <code>{run_dir}</code></p>
<div class="alert">
  <strong>Checkpoint #3:</strong> {branch_note}
  Approve the recommended resolutions or revise via
  <code>Processing-MuAgent revise s7_clustering …</code>.
  See the bottom of this page for exact CLI commands.
</div>

<h2>Recommended resolutions</h2>
<table>
  <tr><th>Modality</th><th>Recommended resolution</th></tr>
  <tr><td>RNA</td><td class="recommended">{rna_res}</td></tr>
  <tr><td>ATAC</td><td class="recommended">{atac_res}</td></tr>
</table>

<h2>RNA sweep</h2>
{rna_table}

<h2>ATAC sweep</h2>
{atac_table}

<h2>Adjacency comparison (tie-breaker)</h2>
{adjacency_block}

<h2>Summary (markdown)</h2>
<pre>{summary_md_escaped}</pre>

<h2>How to approve or revise</h2>
<pre>
# Approve as recommended:
Processing-MuAgent approve s7_clustering --config &lt;run_dir&gt;/deliverables/pre_run/config/run.yaml

# Revise (example — change RNA resolution to 1.2):
Processing-MuAgent revise s7_clustering s7_clustering.rna.resolution=1.2 \\
    --config &lt;run_dir&gt;/deliverables/pre_run/config/run.yaml

# Then resume:
Processing-MuAgent submit --executor pbs   # or 'local' / 'slurm'
</pre>

<p class="meta">
  For interactive exploration (re-run Leiden at different resolutions, view
  side-by-side UMAPs), open <code>resolution_review.ipynb</code> in JupyterLab.
</p>
</body>
</html>
"""


def _sweep_rows_to_html(rows: list[dict[str, Any]], chosen: float | None) -> str:
    if not rows:
        return "<p><em>(sweep not available)</em></p>"
    body = ["<table>",
            "<tr><th>resolution</th><th>n_clusters</th><th>silhouette</th>"
            "<th>stability ARI</th></tr>"]
    for r in rows:
        is_chosen = chosen is not None and abs(float(r["resolution"]) - float(chosen)) < 1e-9
        cls = ' class="recommended"' if is_chosen else ""
        sil = r.get("silhouette", float("nan"))
        body.append(
            f"<tr{cls}>"
            f"<td>{r['resolution']:.2f}</td>"
            f"<td>{r['n_clusters']}</td>"
            f"<td>{sil:.3f}</td>"
            f"<td>{r['seed_stability_ari']:.3f}</td>"
            f"</tr>")
    body.append("</table>")
    return "\n".join(body)


def _adjacency_to_html(report: dict[str, Any]) -> str:
    if not report:
        return "<p><em>(adjacency report not available)</em></p>"
    chunks = []
    for modality, rep in report.items():
        chunks.append(f"<h3>{modality.upper()}</h3>")
        chunks.append(
            f"<p>Recommended: <strong>{rep.get('recommended')}</strong>"
            f"  &nbsp; nearest lower: {rep.get('lower')}"
            f"  &nbsp; nearest higher: {rep.get('higher')}</p>"
        )
        flags = rep.get("surface_flags") or []
        if flags:
            chunks.append('<div class="alert"><strong>Surface flags '
                          '(load-bearing for the decision):</strong><ul>')
            for f in flags:
                chunks.append(f"<li>{html.escape(str(f))}</li>")
            chunks.append("</ul></div>")
        for cmp in rep.get("comparisons", []):
            arrow = "→" if cmp["direction"] == "higher" else "←"
            chunks.append(
                f"<p><code>res={cmp['parent_resolution']}</code> {arrow} "
                f"<code>res={cmp['child_resolution']}</code> &nbsp; "
                f"n_clusters: {cmp['n_parent']} {arrow} {cmp['n_child']} &nbsp; "
                f"ARI: {cmp['ari']} &nbsp; "
                f"<strong>verdict: <code>{cmp['verdict']}</code></strong></p>"
            )
    return "\n".join(chunks)


def build_html(run_dir: Path | str) -> str:
    from ..run_paths import RunPaths
    run_dir = Path(run_dir)
    rp = RunPaths(run_dir)
    art = rp.artifact("s7_clustering")

    import pandas as pd
    import yaml

    sweep = pd.read_parquet(art / "sweep.parquet").to_dict(orient="records")
    rna_rows = [r for r in sweep if r.get("modality") == "rna"]
    atac_rows = [r for r in sweep if r.get("modality") == "atac"]

    params = yaml.safe_load(rp.parameters_yaml.read_text()) or {}
    rna_res = (params.get("s7_clustering.rna.resolution") or {}).get("value")
    atac_res = (params.get("s7_clustering.atac.resolution") or {}).get("value")

    adjacency = {}
    adj_path = art / "adjacency_report.json"
    if adj_path.exists():
        adjacency = json.loads(adj_path.read_text())

    summary_md = ""
    if rp.resolution_summary_md.exists():
        summary_md = rp.resolution_summary_md.read_text()

    return _HTML_TEMPLATE.format(
        run_dir=html.escape(str(run_dir)),
        branch_note=html.escape(_resolution_checkpoint_note(run_dir)),
        rna_res=html.escape(str(rna_res)) if rna_res is not None else "—",
        atac_res=html.escape(str(atac_res)) if atac_res is not None else "—",
        rna_table=_sweep_rows_to_html(rna_rows, rna_res),
        atac_table=_sweep_rows_to_html(atac_rows, atac_res),
        adjacency_block=_adjacency_to_html(adjacency),
        summary_md_escaped=html.escape(summary_md) or "(no summary)",
    )


# ---------------------------------------------------------------------------
# Public entry point — called from s7_clustering_propose rule
# ---------------------------------------------------------------------------

def build_and_render(run_dir: Path | str) -> tuple[Path, Path]:
    """Write resolution_review.{ipynb,html} into deliverables/checkpoint/resolution_review/.

    Returns (ipynb_path, html_path). Safe to call from a local rule on a login
    node — the HTML is rendered statically (no notebook execution required).
    """
    from ..run_paths import RunPaths
    run_dir = Path(run_dir)
    out_dir = RunPaths(run_dir).deliv_resolution_review
    out_dir.mkdir(parents=True, exist_ok=True)

    ipynb = out_dir / "resolution_review.ipynb"
    html_path = out_dir / "resolution_review.html"

    ipynb.write_text(json.dumps(build_notebook(run_dir), indent=1))
    html_path.write_text(build_html(run_dir))

    return ipynb, html_path
