"""Bundle the data that lives next to the package into built wheels.

`config/`, `datasets/schema/` and `openapi.json` sit beside `irv_metrics/` (as in
pawehe's `backend/`). A source checkout or an editable install finds them there;
a built wheel gets a copy under `irv_metrics/_data/`, which `irv_metrics.config`
prefers when present. Everything else is in pyproject.toml.
"""
from pathlib import Path
import shutil

from setuptools import setup
from setuptools.command.build_py import build_py

ROOT = Path(__file__).resolve().parent
DATA = ["config", "datasets/schema", "openapi.json"]


class BuildPyWithData(build_py):
    def run(self):
        super().run()
        if getattr(self, "editable_mode", False):
            return
        target = Path(self.build_lib) / "irv_metrics" / "_data"
        for rel in DATA:
            source, dest = ROOT / rel, target / rel
            if source.is_dir():
                shutil.copytree(source, dest, dirs_exist_ok=True,
                                ignore=shutil.ignore_patterns("__pycache__", "README.md"))
            else:
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, dest)


setup(cmdclass={"build_py": BuildPyWithData})
