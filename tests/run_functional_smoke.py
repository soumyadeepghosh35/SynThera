#!/usr/bin/env python3
"""Real, end-to-end smoke runs of the executables/ scripts against tiny fixture data.

Unlike check_executables_cli.py (which only checks that each script starts up
cleanly), this actually runs DORAnet, DORA-XGB, and ADMET on a handful of
molecules and checks that a real output file comes back. ART is not covered
here - it needs a license and its own training data, well beyond a smoke test
(see build/ART.md); check_executables_cli.py already reports whether it is
installed.

This is slower than the other checks and needs network access the first time
(admet_ai downloads its pretrained model weights on first use). Expect this
to take anywhere from under a minute to several minutes depending on your
machine and connection.

Run from the repository root:

    python tests/run_functional_smoke.py

Environment variables (all optional):
    DORANET_DIR             override the bundled DORAnet checkout path
                             (default: /opt/synthera-deps/doranet)
    DORANET_SMOKE_TIMEOUT   seconds to allow the DORAnet run (default: 600)
    ADMET_SMOKE_TIMEOUT     seconds to allow the ADMET run, including the
                             first-run model download (default: 600)
    DORAXGB_SMOKE_TIMEOUT   seconds to allow the DORA-XGB run (default: 300)

Exit code is 0 if every functional run produced its expected output file.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
EXECUTABLES_DIR = REPO_ROOT / "executables"
FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
OUTPUT_DIR = Path(__file__).resolve().parent / "output"

DORANET_DIR = os.environ.get("DORANET_DIR", "/opt/synthera-deps/doranet")

results = []  # (label, status, detail)


def record(label, status, detail=""):
    results.append((label, status, detail))


def condaRun(env, args, cwd, timeout):
    cmd = ["conda", "run", "-n", env, "python", *args]
    try:
        proc = subprocess.run(
            cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout
        )
        return proc.returncode, proc.stdout, proc.stderr
    except subprocess.TimeoutExpired:
        return None, "", f"timed out after {timeout}s (raise the *_SMOKE_TIMEOUT env var if this needs longer)"


def lastLine(text, fallback="(no output)"):
    stripped = text.strip()
    return stripped.splitlines()[-1] if stripped else fallback


def runDoranetSmoke():
    label = "runDORAnet.py functional run"
    if not Path(DORANET_DIR).is_dir():
        record(label, "FAIL", f"DORAnet checkout not found at {DORANET_DIR}")
        return

    runDir = OUTPUT_DIR / "doranet_smoke"
    runDir.mkdir(parents=True, exist_ok=True)
    jobName = "smoke_ethanol"
    configPath = runDir / "config.yaml"
    configPath.write_text(
        "doranetPath: {doranetPath}\n"
        "jobName: {jobName}\n"
        "starters:\n"
        "  - \"CCO\"\n"
        "generations: 1\n"
        "totalGenerations: 1\n"
        "direction: forward\n"
        "ruleset: JN3604IMT\n"
        "maxAtoms: {{C: 6, N: 3, O: 3, S: 1}}\n".format(doranetPath=DORANET_DIR, jobName=jobName)
    )

    timeout = int(os.environ.get("DORANET_SMOKE_TIMEOUT", "600"))
    returncode, stdout, stderr = condaRun(
        "SynThera_ExptDesigns", [str(EXECUTABLES_DIR / "runDORAnet.py"), str(configPath)],
        cwd=runDir, timeout=timeout,
    )

    outputCsv = runDir / f"{jobName}_molecules.csv"
    if returncode is None:
        record(label, "FAIL", stderr)
    elif returncode != 0:
        record(label, "FAIL", f"exit {returncode} - {lastLine(stdout + stderr)}")
    elif not outputCsv.is_file():
        record(label, "FAIL", f"ran without error but {outputCsv.name} was not created")
    else:
        rowCount = sum(1 for _ in outputCsv.open()) - 1  # minus header
        record(label, "PASS", f"{outputCsv.relative_to(REPO_ROOT)} ({rowCount} molecules)")


def runDoraxgbSmoke():
    label = "runDORAXGB.py functional run"
    fixtureCsv = FIXTURES_DIR / "doraxgb_smoke.csv"
    if not fixtureCsv.is_file():
        record(label, "FAIL", f"fixture not found: {fixtureCsv}")
        return

    runDir = OUTPUT_DIR / "doraxgb_smoke"
    runDir.mkdir(parents=True, exist_ok=True)
    configPath = runDir / "config.yaml"
    configPath.write_text(
        "input:\n"
        "  reactionCSV: \"{reactionCsv}\"\n"
        "output:\n"
        "  fileNamePrefix: \"doraxgb_smoke\"\n"
        "scoring:\n"
        "  testMode: false\n"
        "  useMultiprocessing: false\n"
        "  max_cores: 1\n"
        "  mpChunkSize: 8\n"
        "  cacheSize: 1000\n".format(reactionCsv=fixtureCsv)
    )

    # runDORAXGB.py always writes to reactionDF_wDORAXGBfeasibility.csv in the
    # current working directory, regardless of output.fileNamePrefix - that is
    # a real quirk of the script (the prefix is read but never used), so this
    # test runs from a dedicated directory and checks for that fixed name.
    timeout = int(os.environ.get("DORAXGB_SMOKE_TIMEOUT", "300"))
    returncode, stdout, stderr = condaRun(
        "SynThera_ExptDesigns", [str(EXECUTABLES_DIR / "runDORAXGB.py"), str(configPath)],
        cwd=runDir, timeout=timeout,
    )

    outputCsv = runDir / "reactionDF_wDORAXGBfeasibility.csv"
    if returncode is None:
        record(label, "FAIL", stderr)
    elif returncode != 0:
        record(label, "FAIL", f"exit {returncode} - {lastLine(stdout + stderr)}")
    elif not outputCsv.is_file():
        record(label, "FAIL", f"ran without error but {outputCsv.name} was not created")
    else:
        rowCount = sum(1 for _ in outputCsv.open()) - 1
        record(label, "PASS", f"{outputCsv.relative_to(REPO_ROOT)} ({rowCount} reactions scored)")


def artIsInstalled():
    returncode, _out, _err = condaRun("SynThera", ["-c", "import art.core"], cwd=REPO_ROOT, timeout=30)
    return returncode == 0


def runAdmetSmoke():
    label = "runADMET.py functional run"
    fixtureCsv = FIXTURES_DIR / "admet_smoke.csv"
    if not fixtureCsv.is_file():
        record(label, "FAIL", f"fixture not found: {fixtureCsv}")
        return

    runDir = OUTPUT_DIR / "admet_smoke"
    runDir.mkdir(parents=True, exist_ok=True)
    configPath = runDir / "config.yaml"
    configPath.write_text(
        "inputFile: \"{inputFile}\"\n"
        "outputDir: \"{outputDir}\"\n"
        "SMILES_Column: Canonical_SMILES\n"
        "selectSMILES: null\n"
        "batchSize: 10\n"
        "saveEveryBatch: true\n"
        "resumeFromCheckpoint: false\n"
        "deduplicateSmiles: true\n"
        "enableDiversityFilter: false\n".format(inputFile=fixtureCsv, outputDir=runDir)
    )

    # First run downloads admet_ai's pretrained model weights - needs network
    # and can take a few minutes even for three molecules.
    timeout = int(os.environ.get("ADMET_SMOKE_TIMEOUT", "600"))
    returncode, stdout, stderr = condaRun(
        "SynThera", [str(EXECUTABLES_DIR / "runADMET.py"), str(configPath)],
        cwd=runDir, timeout=timeout,
    )

    outputCsv = runDir / "admet_smoke_wADMET.csv"
    if returncode is None:
        record(label, "FAIL", stderr + " (network access is required for the first run)")
    elif returncode != 0:
        record(label, "FAIL", f"exit {returncode} - {lastLine(stdout + stderr)}")
    elif not outputCsv.is_file():
        record(label, "FAIL", f"ran without error but {outputCsv.name} was not created")
    else:
        rowCount = sum(1 for _ in outputCsv.open()) - 1
        record(label, "PASS", f"{outputCsv.relative_to(REPO_ROOT)} ({rowCount} molecules scored)")


def main():
    if shutil.which("conda") is None:
        print("ERROR: 'conda' is not on PATH. Run this inside the container.")
        sys.exit(2)

    print("\n=== Functional smoke runs (real, tiny data - this can take a while) ===")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Running DORAnet (Step 3) on a single starter molecule...")
    runDoranetSmoke()

    print("Running DORA-XGB (feasibility scoring) on two toy reactions...")
    runDoraxgbSmoke()

    if artIsInstalled():
        print(
            "ART is installed, but a functional smoke test for it is not "
            "included here - it needs training data and a license-specific "
            "setup beyond this suite's scope."
        )

    print("Running ADMET-AI on three molecules (downloads model weights on first run)...")
    runAdmetSmoke()

    labelWidth = max(len(label) for label, _, _ in results) + 2
    print()
    for label, status, detail in results:
        line = f"  [{status:4}] {label:<{labelWidth}}"
        if detail:
            line += f" - {detail}"
        print(line)

    failed = [r for r in results if r[1] == "FAIL"]
    passed = [r for r in results if r[1] == "PASS"]
    print(f"\nfunctional smoke: {len(passed)} passed, {len(failed)} failed (out of {len(results)})")
    print(f"Outputs left under {OUTPUT_DIR.relative_to(REPO_ROOT)}/ for inspection.")

    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
