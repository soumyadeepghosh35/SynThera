# SynThera test suite

This checks that a built SynThera image (or a manual, non-Docker conda setup)
actually has everything the pipeline needs, and that it works - so you don't
find out an environment or a bundled tool is broken partway through a real
Step 1-5 run. It does not validate the science; it validates that the
plumbing is intact: the right packages import in the right environment, and
each `executables/` script actually runs.

## Running it

**If you only have the image** (pulled or built, no source checkout, no
docker-compose bind mount) - the image bakes a standalone copy of
`executables/`, `tests/`, and `external/` into `/opt/synthera-app` for
exactly this case:

```bash
docker run --rm soumya35/synthera:1.0 bash /opt/synthera-app/tests/run_all_tests.sh
docker run --rm soumya35/synthera:1.0 bash /opt/synthera-app/tests/run_all_tests.sh --full
```

Replace `soumya35/synthera:1.0` with whatever tag you built (for example
`synthera:latest` from `docker compose build`). The published image is on
Docker Hub: see the main [README.md](../README.md) for the pull command.
Each invocation is a fresh, throwaway container (`--rm`); output
from a `--full` run is left behind only for the life of that container, so
copy anything you want to keep out with `docker cp` before it exits, or drop
`--rm` and inspect it with `docker cp`/`docker exec` afterward.

**If you have the project checked out** - either running inside the
container via `docker exec -it synthera bash` or a terminal opened inside
JupyterLab (**File -> New -> Terminal**), from the repository root
(`/opt/synthera` in the container), or a running `docker compose` service:

```bash
bash tests/run_all_tests.sh          # fast: tiers 0-2, well under a minute
bash tests/run_all_tests.sh --full   # also tier 3: real runs, several minutes, needs network
```

Both paths run the exact same scripts; the only difference is which copy of
the repository they run against. It prints a `[PASS]` / `[FAIL]` / `[SKIP]`
line per check and exits non-zero if anything required failed. A `[SKIP]` is
not a failure - see "What a SKIP means" below.

## The four tiers

**0. Unit tests** (`unit/`, via `pytest`) - calls individual functions
straight out of `executables/` in-process: config parsing, path building,
validation, and other logic that does not need doranet/rdkit/admet_ai to
actually run anything. No network, no fixture files, done in well under a
second. This is also where a couple of known quirks are pinned down as
regression tests - for example, `runDORAXGB.py`'s `save_output()` silently
ignoring its configured output prefix (see `unit/test_run_doraxgb.py`) - so
a future fix (or regression) shows up here instead of only being noticed
mid pathway design. A module whose own top-level import isn't installed in
the current environment (`DORA_XGB` under `SynThera`, `art.core` unless you
built with `ART_SOURCE`) is skipped there, not failed, via
`pytest.importorskip`.

**1. Environment import checks** (`check_environments.py`) - confirms every
package a notebook or script needs is actually importable in the environment
it is supposed to live in: `SynThera` (Steps 1, 2, 4) and
`SynThera_ExptDesigns` (Steps 3, 5). This is where you'd catch, for example,
a rebuild that silently lost RetroTide, DORAnet, DORA-XGB, or the bundled
`ergochemics` source.

**2. `executables/` CLI smoke checks** (`check_executables_cli.py`) - runs
each of `runDORAnet.py`, `runDORAXGB.py`, `runADMET.py`, and `runART.py` with
no arguments, in its correct conda environment, and confirms it fails
*cleanly* (a usage message) rather than crashing on a missing import. This
takes seconds and catches the two most common breakages: a missing
dependency, or running a script in the wrong environment.

**3. Functional smoke tests** (`run_functional_smoke.py`, opt-in via
`--full`) - actually runs DORAnet, DORA-XGB, and ADMET-AI against a handful of
real molecules from `tests/fixtures/` and checks that a real output file
comes back with the expected number of rows. This is the only tier that
proves the tools genuinely work end to end, not just that they import. It is
slower (DORAnet's runtime depends on how much the ruleset expands even a
single starter molecule; ADMET-AI downloads its pretrained model weights on
first use, which needs network access) - budget several minutes, more on a
slow connection. Outputs are left under `tests/output/` so you can inspect
them afterward.

ART has no functional smoke test: it needs a license and real training data,
which is out of scope for a smoke test (see `build/ART.md`). Tiers 1 and 2
still report whether it is installed.

## What a SKIP means

A `[SKIP]` line is not a failure - it means the check was intentionally
bypassed because the thing it tests is optional in a standard build. The one
you should expect to see on a normal build: **ART is not installed**, in
both `check_environments.py` (env `SynThera`) and `check_executables_cli.py`
(`runART.py`). That's correct - ART is licensed separately and is never
baked into the image by default (see `build/ART.md` for how to add it). If
you did build with `--build-arg ART_SOURCE=...`, this should turn into a real
`[PASS]` instead; if it's still skipped after that, something about the ART
install didn't take.

## If something actually fails

- **A failure in tier 0** (pytest reports it directly, not as a `[FAIL]`
  line) means the logic inside an `executables/` script itself is broken -
  a bad edit to config parsing, path building, or similar. Read pytest's own
  assertion output; it names the exact function and expectation that broke.
- **A `[FAIL]` in tier 1** for a specific package means that package isn't
  importable in the environment that's supposed to provide it. Check the
  detail message (usually `ModuleNotFoundError`) and cross-reference
  `build/environment.yml` / `build/exptDesigns_environment.yml` and
  `build/Dockerfile` - either the build didn't finish cleanly, or the image
  you're running predates a fix. Rebuilding with `docker compose up --build`
  is the first thing to try.
- **A `[FAIL]` in tier 2** with "unexpected traceback" means the script
  crashed on start rather than showing its usage message - almost always the
  same root cause as a tier 1 failure for a package that script depends on.
- **A `[FAIL]` in tier 3** means the script ran but didn't produce the output
  it should have, or crashed partway through. The detail line includes the
  script's own last line of output, which usually says why (for ADMET, a
  network error trying to download model weights is the most common cause).

## Files here

- `run_all_tests.sh` - the orchestrator; run this, not the individual scripts,
  unless you're debugging one check in isolation.
- `unit/` - tier 0, run with `pytest` once per environment (the shell script
  already does this for both); see `unit/conftest.py` for how it locates
  `executables/`.
- `check_environments.py` - tier 1, run once per environment (the shell
  script already does this for both).
- `check_executables_cli.py` - tier 2, run once (it shells out to both
  environments itself).
- `run_functional_smoke.py` - tier 3, run once (same).
- `fixtures/` - the tiny, deliberately trivial input files tier 3 uses
  (three drug-like molecules for ADMET, two toy reactions for DORA-XGB;
  DORAnet's tier 3 test generates its own one-line config against ethanol).
- `output/` - created on first `--full` run; safe to delete between runs.
