"""Reproduce the SynThera conda environments and register their Jupyter kernels.

This script builds both environments required by the SynThera pipeline from their
exported specification files and registers a matching Jupyter kernel for each, so
that a user can select the correct interpreter for every notebook.

  environment.yml            -> environment 'SynThera'             (Steps 1-4)
  exptDesigns_environment.yml -> environment 'SynThera_ExptDesigns' (Step 5)

Usage:
  python setup.py                 # create both environments and kernels
  python setup.py --force         # remove and recreate environments that exist
  python setup.py --skip-kernels  # create environments only
  python setup.py --only SynThera # act on a single target
"""

import argparse
import platform
import shutil
import subprocess
import sys
from pathlib import Path

CONFIG = [
    {
        "ymlFile": "environment.yml",
        "envName": "SynThera",
        "kernelDisplayName": "SynThera",
    },
    {
        "ymlFile": "exptDesigns_environment.yml",
        "envName": "SynThera_ExptDesigns",
        "kernelDisplayName": "SynThera_ExptDesigns",
    },
]

SCRIPT_DIR = Path(__file__).resolve().parent


def parseArgs():
    parser = argparse.ArgumentParser(
        description="Build the SynThera conda environments and Jupyter kernels."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Remove and recreate an environment if it already exists.",
    )
    parser.add_argument(
        "--skip-kernels",
        dest="skipKernels",
        action="store_true",
        help="Create the environments without registering Jupyter kernels.",
    )
    parser.add_argument(
        "--only",
        choices=[entry["envName"] for entry in CONFIG],
        help="Restrict the run to a single environment.",
    )
    return parser.parse_args()


def findCondaExecutable():
    """Return the preferred solver executable, favouring mamba for speed."""
    for candidate in ("mamba", "conda"):
        if shutil.which(candidate):
            return candidate
    raise RuntimeError(
        "Neither 'mamba' nor 'conda' was found on PATH. Install Miniconda or "
        "Miniforge first, then rerun this script."
    )


def runCommand(command, check=True):
    """Run a command, streaming its output to the terminal."""
    print(f"\n$ {' '.join(command)}\n")
    result = subprocess.run(command, check=False)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}: {' '.join(command)}"
        )
    return result.returncode


def checkPlatform(condaExe):
    """Warn if the host is not Linux, since the specs are build-pinned linux-64."""
    if platform.system() != "Linux":
        print(
            "\nWarning: these environment files are fully build-pinned linux-64 "
            "exports. On "
            f"{platform.system()} the exact builds will not solve. Use a Linux "
            "machine or container, or recreate the specs from history on your "
            "platform. See the README for details.\n"
        )


def environmentExists(condaExe, envName):
    result = subprocess.run(
        [condaExe, "env", "list"],
        check=True,
        capture_output=True,
        text=True,
    )
    for line in result.stdout.splitlines():
        token = line.split()
        if token and token[0] == envName:
            return True
    return False


def createEnvironment(condaExe, spec, force):
    envName = spec["envName"]
    ymlPath = SCRIPT_DIR / spec["ymlFile"]

    if not ymlPath.is_file():
        print(
            f"Skipping '{envName}': specification file '{ymlPath.name}' was not "
            f"found in {SCRIPT_DIR}. Place it beside this script and rerun."
        )
        return False

    if environmentExists(condaExe, envName):
        if force:
            print(f"Environment '{envName}' exists. Removing it because --force was set.")
            runCommand([condaExe, "env", "remove", "-n", envName, "-y"])
        else:
            print(
                f"Environment '{envName}' already exists. Leaving it untouched. "
                "Pass --force to recreate it."
            )
            return True

    print(f"Creating environment '{envName}' from '{ymlPath.name}'.")
    runCommand([condaExe, "env", "create", "-n", envName, "-f", str(ymlPath)])
    return True


def registerKernel(condaExe, spec):
    """Install ipykernel if needed, then register a kernel named after the env."""
    envName = spec["envName"]
    displayName = spec["kernelDisplayName"]

    hasKernel = subprocess.run(
        [condaExe, "run", "-n", envName, "python", "-c", "import ipykernel"],
        check=False,
        capture_output=True,
    )
    if hasKernel.returncode != 0:
        print(f"ipykernel is missing in '{envName}'. Installing it now.")
        runCommand([condaExe, "install", "-n", envName, "-y", "ipykernel"])

    print(f"Registering Jupyter kernel '{envName}'.")
    runCommand(
        [
            condaExe,
            "run",
            "-n",
            envName,
            "python",
            "-m",
            "ipykernel",
            "install",
            "--user",
            "--name",
            envName,
            "--display-name",
            displayName,
        ]
    )


def main():
    args = parseArgs()
    condaExe = findCondaExecutable()
    print(f"Using solver: {condaExe}")
    checkPlatform(condaExe)

    targets = CONFIG
    if args.only:
        targets = [entry for entry in CONFIG if entry["envName"] == args.only]

    created = []
    for spec in targets:
        wasCreated = createEnvironment(condaExe, spec, args.force)
        if wasCreated:
            created.append(spec)

    if not args.skipKernels:
        for spec in created:
            registerKernel(condaExe, spec)

    print("\nDone.")
    if created:
        names = ", ".join(spec["envName"] for spec in created)
        print(f"Environments ready: {names}")
    if not args.skipKernels and created:
        print("Kernels are visible in Jupyter under the same names.")
    print("Activate an environment with: conda activate <name>")


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as error:
        print(f"\nError: {error}", file=sys.stderr)
        sys.exit(1)
