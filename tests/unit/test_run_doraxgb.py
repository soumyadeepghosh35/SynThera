"""Unit tests for the pure-logic pieces of executables/runDORAXGB.py.

runDORAXGB.py imports DORA_XGB at module level, and that package is only
installed into SynThera_ExptDesigns (see build/Dockerfile) - so this whole
file is skipped when run under SynThera.
"""

import pandas as pd
import pytest

runDORAXGB = pytest.importorskip(
    "runDORAXGB",
    reason="DORA-XGB is only installed in SynThera_ExptDesigns, not SynThera",
)


def testLoadConfigReadsYaml(tmp_path):
    configPath = tmp_path / "config.yaml"
    configPath.write_text("scoring:\n  max_cores: 2\n")
    assert runDORAXGB.load_config(configPath) == {"scoring": {"max_cores": 2}}


def testLoadReactionsRequiresReactionStringColumn(tmp_path):
    csvPath = tmp_path / "reactions.csv"
    csvPath.write_text("notTheRightColumn\nfoo\n")
    with pytest.raises(ValueError, match="reactionString"):
        runDORAXGB.load_reactions(csvPath)


def testLoadReactionsStripsSpacesAndRespectsTestMode(tmp_path):
    csvPath = tmp_path / "reactions.csv"
    csvPath.write_text("reactionString\nCC O>>CC = O\nCC(=O)O>>CC(=O)OCC\n")
    df = runDORAXGB.load_reactions(csvPath, test_mode=True, test_rows=1)
    assert len(df) == 1
    assert df.iloc[0]["rxn_str"] == "CCO>>CC=O"


def testGetDefaultCoreIdsRaisesWhenTooManyRequested():
    with pytest.raises(ValueError, match="max_cores"):
        runDORAXGB.get_default_core_ids(max_cores=10_000)


def testSaveOutputIgnoresPrefixAndWritesFixedFilename(tmp_path, monkeypatch):
    """Pins a known quirk: save_output takes a prefix argument but never uses
    it - the output always lands at a fixed filename in the current working
    directory (see build/DOCKER.md and tests/run_functional_smoke.py, which
    both work around this). If this is ever fixed upstream, update this test
    and those docs to match the new, prefix-aware behavior.
    """
    monkeypatch.chdir(tmp_path)
    runDORAXGB.save_output(
        pd.DataFrame({"rxn_str": ["CCO>>CC=O"]}), prefix="should_be_ignored"
    )

    fixedOutput = tmp_path / "reactionDF_wDORAXGBfeasibility.csv"
    prefixedOutput = tmp_path / "should_be_ignored.csv"
    assert fixedOutput.is_file()
    assert not prefixedOutput.is_file()
