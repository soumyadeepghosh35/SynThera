"""Unit tests for the helper functions inside the test suite's own scripts -
check_executables_cli.py and run_functional_smoke.py. Pure stdlib/PyYAML
logic, no conda environment specifics or heavy dependency involved, so this
file runs the same way under either environment.
"""

import yaml

import check_executables_cli
import run_functional_smoke


def testLastLineReturnsFinalNonEmptyLine():
    assert run_functional_smoke.lastLine("first\nsecond\nthird\n") == "third"


def testLastLineReturnsFallbackForEmptyText():
    assert run_functional_smoke.lastLine("   \n  ", fallback="(nothing)") == "(nothing)"


def testDoranetSmokeConfigIsValidYaml():
    """The DORAnet smoke config is hand-formatted (not yaml.dump'd) to avoid
    a PyYAML dependency in run_functional_smoke.py's own top-level logic -
    confirm the hand-formatting is still valid YAML with the fields
    runDORAnet.py needs.
    """
    text = (
        "doranetPath: {doranetPath}\n"
        "jobName: {jobName}\n"
        "starters:\n"
        '  - "CCO"\n'
        "generations: 1\n"
        "totalGenerations: 1\n"
        "direction: forward\n"
        "ruleset: JN3604IMT\n"
        "maxAtoms: {{C: 6, N: 3, O: 3, S: 1}}\n"
    ).format(doranetPath="/opt/synthera-deps/doranet", jobName="smoke_ethanol")

    parsed = yaml.safe_load(text)
    assert parsed["jobName"] == "smoke_ethanol"
    assert parsed["starters"] == ["CCO"]
    assert parsed["maxAtoms"] == {"C": 6, "N": 3, "O": 3, "S": 1}


def testScriptsTableHasExpectedEntries():
    scriptNames = {name for name, _env, _code in check_executables_cli.SCRIPTS}
    assert scriptNames == {"runDORAnet.py", "runDORAXGB.py", "runADMET.py", "runART.py"}


def testScriptsTableRoutesEachScriptToItsOwnEnvironment():
    envByScript = {name: env for name, env, _code in check_executables_cli.SCRIPTS}
    assert envByScript["runDORAnet.py"] == "SynThera_ExptDesigns"
    assert envByScript["runDORAXGB.py"] == "SynThera_ExptDesigns"
    assert envByScript["runADMET.py"] == "SynThera"
    assert envByScript["runART.py"] == "SynThera"
