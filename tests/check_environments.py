#!/usr/bin/env python3
"""Import-level health check for one SynThera conda environment.

Confirms that the packages a given environment is supposed to provide are
actually importable - including the three tools baked in from their own
GitHub sources (RetroTide, DORAnet, DORA-XGB) and the bundled ergochemics
library, which are not ordinary pip installs. This does not run any real
compute; it only checks that every import SynThera's notebooks and
executables/ scripts rely on actually resolves in this environment.

Run inside each target environment (not the base environment):

    conda run -n SynThera             python tests/check_environments.py SynThera
    conda run -n SynThera_ExptDesigns python tests/check_environments.py SynThera_ExptDesigns

Exit code is 0 if every required check passed (a SKIP, such as ART not being
installed, does not count as a failure); non-zero if anything required failed.
"""

import importlib
import os
import sys
from pathlib import Path


# tests/ lives one level below the repository root in both places this file
# ever runs from: the live, bind-mounted repo at /opt/synthera, and the
# standalone copy baked into the image at /opt/synthera-app for a user who
# only has the image (no full source checkout, so no Input_Data directory to
# search for). Walking up from this file's own location works identically in
# both cases and does not depend on any directory that a slimmed-down copy
# might not have.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DORANET_DIR = os.environ.get("DORANET_DIR", "/opt/synthera-deps/doranet")
ERGOCHEMICS_SRC = os.environ.get(
    "ERGOCHEMICS_SRC", str(PROJECT_ROOT / "external" / "ergochemics" / "src")
)

results = []  # list of (label, status, detail)


def record(label, ok, detail="", skip=False):
    status = "SKIP" if skip else ("PASS" if ok else "FAIL")
    results.append((label, status, detail))


def checkPlainImport(moduleName, label=None, requiredIfMissing=True):
    label = label or moduleName
    try:
        importlib.import_module(moduleName)
        record(label, True)
    except Exception as exc:  # noqa: BLE001 - we want to report any failure
        record(label, False, f"{type(exc).__name__}: {exc}", skip=not requiredIfMissing)


def checkDoranet():
    label = f"doranet (bundled checkout at {DORANET_DIR})"
    try:
        if not Path(DORANET_DIR).is_dir():
            raise FileNotFoundError(f"DORAnet checkout not found at {DORANET_DIR}")
        if DORANET_DIR not in sys.path:
            sys.path.insert(0, DORANET_DIR)
        importlib.import_module("doranet.modules.enzymatic")
        importlib.import_module("doranet.modules.post_processing")
        record(label, True)
    except Exception as exc:  # noqa: BLE001
        record(label, False, f"{type(exc).__name__}: {exc}")


def checkErgochemics():
    label = f"ergochemics (bundled source at {ERGOCHEMICS_SRC})"
    try:
        if not Path(ERGOCHEMICS_SRC).is_dir():
            raise FileNotFoundError(f"ergochemics source not found at {ERGOCHEMICS_SRC}")
        if ERGOCHEMICS_SRC not in sys.path:
            sys.path.insert(0, ERGOCHEMICS_SRC)
        importlib.import_module("ergochemics.mapping")
        importlib.import_module("ergochemics.similarity")
        record(label, True)
    except Exception as exc:  # noqa: BLE001
        record(label, False, f"{type(exc).__name__}: {exc}")


def checkRetrotide():
    try:
        importlib.import_module("retrotide")
        importlib.import_module("bcs")
        record("retrotide + bcs (RetroTide)", True)
    except Exception as exc:  # noqa: BLE001
        record("retrotide + bcs (RetroTide)", False, f"{type(exc).__name__}: {exc}")


def checkDoraXgb():
    # Same import runExecutables/runDORAXGB.py uses: `from DORA_XGB import DORA_XGB`
    # gives the DORA_XGB submodule, which is what provides feasibility_classifier().
    try:
        from DORA_XGB import DORA_XGB as doraXgbModule
        if not callable(getattr(doraXgbModule, "feasibility_classifier", None)):
            raise ImportError("DORA_XGB.DORA_XGB has no feasibility_classifier()")
        record("DORA_XGB (DORA-XGB)", True)
    except Exception as exc:  # noqa: BLE001
        record("DORA_XGB (DORA-XGB)", False, f"{type(exc).__name__}: {exc}")


def checkArt():
    """ART is optional and licensed separately - a missing import is a SKIP, not a FAIL."""
    try:
        importlib.import_module("art.core")
        record("art (ART - optional, licensed)", True)
    except Exception as exc:  # noqa: BLE001
        record(
            "art (ART - optional, licensed)",
            False,
            f"not installed ({type(exc).__name__}) - expected unless built with "
            "--build-arg ART_SOURCE=... ; see build/ART.md",
            skip=True,
        )


SYNTHERA_CHECKS = [
    lambda: checkPlainImport("rdkit"),
    lambda: checkPlainImport("pandas"),
    lambda: checkPlainImport("numpy"),
    lambda: checkPlainImport("scipy"),
    lambda: checkPlainImport("sklearn", label="scikit-learn"),
    lambda: checkPlainImport("torch"),
    lambda: checkPlainImport("tensorflow"),
    lambda: checkPlainImport("transformers"),
    lambda: checkPlainImport("chemprop"),
    lambda: checkPlainImport("xgboost"),
    lambda: checkPlainImport("macaw"),
    lambda: checkPlainImport("admet_ai"),
    checkDoranet,
    checkArt,
]

SYNTHERA_EXPTDESIGNS_CHECKS = [
    lambda: checkPlainImport("rdkit"),
    lambda: checkPlainImport("pandas"),
    lambda: checkPlainImport("numpy"),
    lambda: checkPlainImport("cobra"),
    lambda: checkPlainImport("dnachisel"),
    lambda: checkPlainImport("Bio", label="Bio (biopython)"),
    lambda: checkPlainImport("equilibrator_api"),
    lambda: checkPlainImport("xgboost"),
    lambda: checkPlainImport("sklearn", label="scikit-learn"),
    lambda: checkPlainImport("mapchiral"),
    checkRetrotide,
    checkDoranet,
    checkDoraXgb,
    checkErgochemics,
]

ENV_CHECKS = {
    "SynThera": SYNTHERA_CHECKS,
    "SynThera_ExptDesigns": SYNTHERA_EXPTDESIGNS_CHECKS,
}


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in ENV_CHECKS:
        print(f"Usage: python {Path(__file__).name} <SynThera|SynThera_ExptDesigns>")
        sys.exit(2)

    envName = sys.argv[1]
    actualEnv = os.environ.get("CONDA_DEFAULT_ENV", "<unknown - not run via conda run>")

    print(f"\n=== Environment check: {envName} ===")
    print(f"(conda reports the active environment as: {actualEnv})")
    if actualEnv != envName:
        print(
            f"WARNING: you asked to check '{envName}' but the active conda "
            f"environment is '{actualEnv}'. Run this with "
            f"`conda run -n {envName} python {Path(__file__).name} {envName}`."
        )

    for check in ENV_CHECKS[envName]:
        check()

    labelWidth = max(len(label) for label, _, _ in results) + 2
    for label, status, detail in results:
        line = f"  [{status:4}] {label:<{labelWidth}}"
        if detail and status != "PASS":
            line += f" - {detail}"
        print(line)

    failed = [r for r in results if r[1] == "FAIL"]
    skipped = [r for r in results if r[1] == "SKIP"]
    passed = [r for r in results if r[1] == "PASS"]
    print(
        f"\n{envName}: {len(passed)} passed, {len(skipped)} skipped, "
        f"{len(failed)} failed (out of {len(results)})"
    )

    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
