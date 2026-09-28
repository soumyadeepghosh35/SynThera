"""Shared setup for the SynThera unit test tier (tests/unit/).

Unlike check_environments.py (import-level health) and
check_executables_cli.py / run_functional_smoke.py (run each script as a
real subprocess), this tier calls individual functions in-process: config
parsing, path building, validation, and similar logic that does not need a
live compute run. That makes it fast and network-free enough to run on
every build.

Run inside a conda environment that provides whatever the module under test
imports at its own top level - see tests/README.md for which environment
each test file expects:

    conda run -n SynThera             pytest tests/unit
    conda run -n SynThera_ExptDesigns pytest tests/unit

A module whose own top-level import needs a tool that is not installed in
the current environment (DORA_XGB, art.core) is skipped there via
pytest.importorskip, not failed - the same PASS/FAIL/SKIP convention used
throughout this test suite.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
EXECUTABLES_DIR = REPO_ROOT / "executables"
TESTS_DIR = REPO_ROOT / "tests"

for directory in (EXECUTABLES_DIR, TESTS_DIR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
