#!/usr/bin/env bash
# SynThera test suite - confirms the Docker image (or a manual conda setup)
# actually has everything the pipeline needs, working, before you trust it
# for a real run.
#
# Works two ways:
#   - With a source checkout / docker-compose bind mount, from the repo root
#     (/opt/synthera inside the container):
#       bash tests/run_all_tests.sh
#       bash tests/run_all_tests.sh --full
#   - From just the built/pulled image, no bind mount at all:
#       docker run --rm <image> bash /opt/synthera-app/tests/run_all_tests.sh
#       docker run --rm <image> bash /opt/synthera-app/tests/run_all_tests.sh --full
#   - Or against a container that is already running:
#       docker exec -it synthera bash tests/run_all_tests.sh
#       docker exec -it synthera bash tests/run_all_tests.sh --full
#
# Exit code is 0 only if nothing failed. A SKIP (most commonly: ART is not
# installed, which is expected unless you built with --build-arg
# ART_SOURCE=...) is not treated as a failure. See tests/README.md.

set -uo pipefail

cd "$(dirname "$0")/.." || exit 1
REPO_ROOT="$(pwd)"

echo "=================================================================="
echo " SynThera test suite"
echo " repo root: $REPO_ROOT"
echo "=================================================================="

overallStatus=0

echo
echo "---- 0. Unit tests (function-level, no compute run) ----------------"
conda run -n SynThera             pytest -q tests/unit || overallStatus=1
conda run -n SynThera_ExptDesigns pytest -q tests/unit || overallStatus=1

echo
echo "---- 1. Environment import checks --------------------------------"
conda run -n SynThera             python tests/check_environments.py SynThera             || overallStatus=1
conda run -n SynThera_ExptDesigns python tests/check_environments.py SynThera_ExptDesigns || overallStatus=1

echo
echo "---- 2. executables/ CLI smoke checks ------------------------------"
python tests/check_executables_cli.py || overallStatus=1

if [ "${1:-}" = "--full" ]; then
  echo
  echo "---- 3. Functional smoke tests (slow, needs network) ---------------"
  python tests/run_functional_smoke.py || overallStatus=1
else
  echo
  echo "(skipping functional smoke tests - rerun with --full to include real,"
  echo " tiny end-to-end runs of DORAnet, DORA-XGB, and ADMET)"
fi

echo
echo "=================================================================="
if [ "$overallStatus" -eq 0 ]; then
  echo " RESULT: PASS - see any SKIP lines above for what was intentionally left out"
else
  echo " RESULT: FAIL - see the [FAIL] lines above for what to fix"
fi
echo "=================================================================="

exit "$overallStatus"
