"""Regression test: the plan-review intro paragraph must survive re-renders.

The agent writes the intro via `executor plan-review --intro "<prose>"`, but the
`plan_review_propose` Snakemake rule re-renders plan_review.md AND
plan_summary.html with no intro argument. The intro must be persisted so those
intro-less re-renders reproduce it instead of silently dropping it.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from executor import plan_review as pr
from executor.run_paths import RunPaths

_INTRO = "UNIQUE_INTRO_SENTINEL describing a paired multiome testis sample."


class PlanIntroPersistTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.run_dir = Path(self._tmp.name)
        self.paths = RunPaths(self.run_dir)
        self.paths.ensure()
        # Pre-create an (empty) qc_explore.json so render does not try to compute it.
        qc = self.paths.artifact("qc_explore", "qc_explore.json")
        qc.parent.mkdir(parents=True, exist_ok=True)
        qc.write_text("{}")

    def tearDown(self):
        self._tmp.cleanup()

    def test_helper_roundtrip(self):
        pr._persist_intro(self.paths, _INTRO)
        self.assertEqual(pr._load_persisted_intro(self.paths), _INTRO)
        self.assertTrue(self.paths.plan_intro.exists())

    def test_intro_survives_introless_rerender(self):
        # 1. Agent renders with --intro: persists + writes both deliverables.
        pr.write_summary(self.run_dir, intro=_INTRO)
        pr.write_plan_summary_html(self.run_dir, intro=_INTRO)
        self.assertIn(_INTRO, self.paths.plan_review_md.read_text())
        self.assertIn(_INTRO, self.paths.plan_summary_html.read_text())

        # 2. Propose rule re-renders with NO intro (simulates plan_review_propose).
        pr.write_summary(self.run_dir)
        pr.write_plan_summary_html(self.run_dir)

        # 3. Intro must still be present in BOTH deliverables.
        self.assertIn(_INTRO, self.paths.plan_review_md.read_text(),
                      "intro dropped from plan_review.md on intro-less re-render")
        self.assertIn(_INTRO, self.paths.plan_summary_html.read_text(),
                      "intro dropped from plan_summary.html on intro-less re-render")

    def test_no_persisted_intro_renders_without_error(self):
        # No intro ever set: render must still succeed (and contain no sentinel).
        pr.write_summary(self.run_dir)
        self.assertNotIn(_INTRO, self.paths.plan_review_md.read_text())


if __name__ == "__main__":
    unittest.main()
