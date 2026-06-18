"""launch_runner.sh must put the repo on PYTHONPATH so every job tier imports `executor`
from source (parity with the GPU container wrapper + the head job's cwd=repo_root). This
makes a submit-time auto-provisioned env (created from the lock, no `pip install -e`)
sufficient for CPU child jobs, which inherit PYTHONPATH via sbatch/qsub --export=ALL.
"""
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LAUNCH = REPO / "scripts" / "launch_runner.sh"


class LaunchRunnerPythonPathTests(unittest.TestCase):
    def test_exports_pythonpath_with_repo_root(self):
        lines = LAUNCH.read_text().splitlines()
        self.assertTrue(
            any("export PYTHONPATH=" in ln and "PMA_REPO_ROOT" in ln for ln in lines),
            "launch_runner.sh must export PYTHONPATH including $PMA_REPO_ROOT so child "
            "jobs import executor without requiring pip install -e")


if __name__ == "__main__":
    unittest.main()
