"""Filter the exported SynThera environment specifications for container builds.

The exported specs contain pip entries that either cannot resolve from a clean
PyPI build, or that the Dockerfile deliberately reinstalls from their own
GitHub sources afterward so the exact code the executables/ scripts expect is
guaranteed:

  art              Automated Recommendation Tool. Excluded on purpose: ART is
                   licensed separately for academic use, so it is never baked
                   into the image. Install it yourself after obtaining a license
                   and the author's instructions (see DOCKER.md).
  pathermo         Not published to PyPI. Removed so the build can proceed;
                   supply your own source if a downstream step needs it.
  doranet          Reinstalled by the Dockerfile from
                   https://github.com/wsprague-nu/doranet.git into both
                   environments, from a checkout kept on disk in the image
                   (runDORAnet.py needs a real directory, not just an
                   installed package). Removed here so the PyPI pin is never
                   installed first only to be overwritten.
  dora-xgb         Reinstalled by the Dockerfile from
                   https://github.com/tyo-nu/DORA_XGB.git into
                   SynThera_ExptDesigns only, for the same reason.
  python-graphviz  The conda name for the binding published on PyPI as
                   'graphviz'. Renamed so pip can resolve it.

The filter edits each spec line by line, so every build-pinned conda and pip
entry is preserved byte for byte. Run it before 'conda env create'.
"""

import sys
from pathlib import Path

FILTERS = {
    "environment.yml": {
        "remove": ("art==", "pathermo==", "doranet=="),
        "rename": {},
    },
    "exptDesigns_environment.yml": {
        "remove": ("pathermo==", "doranet==", "dora-xgb=="),
        "rename": {"python-graphviz==0.21": "graphviz==0.21"},
    },
}


def pipToken(line):
    """Return the dependency token from a pip list line, or None."""
    stripped = line.strip()
    if stripped.startswith("- "):
        return stripped[2:].strip()
    return None


def filterSpec(path, removePrefixes, renameMap):
    """Rewrite one spec file in place, returning counts of the edits made."""
    original = path.read_text().splitlines(keepends=True)
    kept = []
    removed = 0
    renamed = 0

    for line in original:
        token = pipToken(line)
        if token is not None:
            if any(token.startswith(prefix) for prefix in removePrefixes):
                print(f"  remove: {token}")
                removed += 1
                continue
            if token in renameMap:
                replacement = renameMap[token]
                print(f"  rename: {token} -> {replacement}")
                line = line.replace(token, replacement, 1)
                renamed += 1
        kept.append(line)

    path.write_text("".join(kept))
    return removed, renamed


def main():
    scriptDir = Path(__file__).resolve().parent
    totalRemoved = 0
    totalRenamed = 0

    for fileName, rules in FILTERS.items():
        path = scriptDir / fileName
        if not path.is_file():
            print(f"Skipping '{fileName}': not found in {scriptDir}.")
            continue
        print(f"Filtering {fileName}")
        removed, renamed = filterSpec(path, rules["remove"], rules["rename"])
        totalRemoved += removed
        totalRenamed += renamed

    print(f"\nDone. Removed {totalRemoved} entrie(s), renamed {totalRenamed}.")
    print("ART is intentionally excluded; install it separately once licensed.")
    print("DORAnet and DORA-XGB are reinstalled by the Dockerfile from their own")
    print("GitHub sources; see this file's module docstring.")


if __name__ == "__main__":
    try:
        main()
    except OSError as error:
        print(f"\nError: {error}", file=sys.stderr)
        sys.exit(1)
