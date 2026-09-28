#!/usr/bin/env python3
"""Fast CLI smoke test for every script in executables/.

This does not run any real compute - it only confirms that each script's own
module-level imports succeed in its correct conda environment, and that it
fails *cleanly* (a usage message, not a traceback) when called without a
config file. That is enough to catch the two most common ways an executables/
script breaks after a rebuild: a missing dependency, or a wrong conda
environment.

For real, end-to-end runs against tiny sample data, see run_functional_smoke.py
(slower, needs network for the ADMET model download).

Run from the repository root, from any environment that has `conda` on PATH
(the base environment JupyterLab starts in is fine):

    python tests/check_executables_cli.py

Exit code is 0 if every script behaved as expected (ART being unavailable is
reported as SKIP, not FAIL, since it is optional and licensed separately).
"""

import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
EXECUTABLES_DIR = REPO_ROOT / "executables"
TIMEOUT_SECONDS = 90

# (script, conda environment, expected "clean failure" exit code when run
#  with no arguments). runDORAnet.py prints its own usage message and exits
#  1; the other three use argparse, whose own usage error exits 2.
SCRIPTS = [
    ("runDORAnet.py", "SynThera_ExptDesigns", 1),
    ("runDORAXGB.py", "SynThera_ExptDesigns", 2),
    ("runADMET.py", "SynThera", 2),
    ("runART.py", "SynThera", 2),
]

results = []  # (label, status, detail)


def record(label, status, detail=""):
    results.append((label, status, detail))


def condaRun(env, args, timeout=TIMEOUT_SECONDS):
    cmd = ["conda", "run", "-n", env, "python", *args]
    try:
        proc = subprocess.run(
            cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=timeout
        )
        return proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired:
        return None, "", f"timed out after {timeout}s"


def artIsInstalled():
    returncode, _out, _err = condaRun("SynThera", ["-c", "import art.core"], timeout=30)
    return returncode == 0


def checkScript(scriptName, env, expectedExitCode):
    scriptPath = EXECUTABLES_DIR / scriptName
    label = f"{scriptName} ({env})"

    if not scriptPath.is_file():
        record(label, "FAIL", f"file not found: {scriptPath}")
        return

    if scriptName == "runART.py" and not artIsInstalled():
        record(
            label,
            "SKIP",
            "ART not installed - expected unless built with "
            "--build-arg ART_SOURCE=... ; see build/ART.md",
        )
        return

    returncode, stdout, stderr = condaRun(env, [str(scriptPath)])

    if returncode is None:
        record(label, "FAIL", stderr)
        return

    combined = stdout + stderr
    hasTraceback = "Traceback (most recent call last)" in combined

    if hasTraceback:
        firstErrorLine = next(
            (line for line in combined.splitlines() if line.strip().startswith(
                ("ModuleNotFoundError", "ImportError", "FileNotFoundError")
            )),
            combined.strip().splitlines()[-1] if combined.strip() else "(no output)",
        )
        record(label, "FAIL", f"unexpected traceback - {firstErrorLine}")
    elif returncode == expectedExitCode:
        record(label, "PASS", f"exited {returncode} with a clean usage message, as expected")
    else:
        detail = (combined.strip().splitlines() or ["(no output)"])[-1]
        record(
            label,
            "FAIL",
            f"expected exit code {expectedExitCode}, got {returncode} - {detail}",
        )


def main():
    if shutil.which("conda") is None:
        print(
            "ERROR: 'conda' is not on PATH. Run this from a shell where conda "
            "is available (for example, a JupyterLab terminal in the base "
            "environment, or `docker exec -it synthera bash`)."
        )
        sys.exit(2)

    print("\n=== executables/ CLI smoke test ===")
    for scriptName, env, expectedExitCode in SCRIPTS:
        checkScript(scriptName, env, expectedExitCode)

    labelWidth = max(len(label) for label, _, _ in results) + 2
    for label, status, detail in results:
        line = f"  [{status:4}] {label:<{labelWidth}}"
        if detail:
            line += f" - {detail}"
        print(line)

    failed = [r for r in results if r[1] == "FAIL"]
    skipped = [r for r in results if r[1] == "SKIP"]
    passed = [r for r in results if r[1] == "PASS"]
    print(
        f"\nexecutables/: {len(passed)} passed, {len(skipped)} skipped, "
        f"{len(failed)} failed (out of {len(results)})"
    )

    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
