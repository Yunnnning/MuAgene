"""`regenerate-locks` must reject a `pip:` subsection in the CPU env YAML.

The CPU lock is rendered with `conda-lock --kind explicit` (conda-only). A `pip:` block
would be silently dropped from the lock, so a freshly provisioned env would be missing
those packages and fail `validate-env` far from this command. The guard turns that silent
data loss into a loud, local error. A bare `- pip` string (needed for editable agent
installs) must still be allowed.
"""
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import yaml
from click.testing import CliRunner

from executor import cli


def _manifest():
    return {"cpu": {"definition": "env.yaml", "lock": "env.lock"}}


class RegenerateLocksPipGuardTests(unittest.TestCase):
    def _run(self, yaml_text):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "env.yaml").write_text(yaml_text)
            with mock.patch.object(cli.hpc, "REPO_ROOT", Path(tmp)), \
                 mock.patch.object(cli.hpc, "load_env_manifest", return_value=_manifest()):
                return CliRunner().invoke(cli.main, ["regenerate-locks"])

    def test_pip_subsection_fails_loud(self):
        res = self._run(
            "name: muagene\nchannels: [conda-forge]\n"
            "dependencies:\n  - python\n  - pip\n  - pip:\n      - scrublet\n")
        self.assertNotEqual(res.exit_code, 0)
        self.assertIn("pip:", res.output)
        self.assertIn("silently dropped", res.output)

    def test_bare_pip_string_passes_guard(self):
        # No pip: mapping -> guard passes; the command then stops at the conda-lock check
        # (conda-lock may be absent in CI). Either way it must NOT hit the pip guard.
        res = self._run(
            "name: muagene\nchannels: [conda-forge]\ndependencies:\n  - python\n  - pip\n")
        self.assertNotIn("silently dropped", res.output)


if __name__ == "__main__":
    unittest.main()
