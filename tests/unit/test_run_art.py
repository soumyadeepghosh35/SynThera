"""Unit tests for the pure-logic pieces of executables/runART.py.

runART.py imports art.core at module level, and ART is licensed separately
and not installed by default (see build/ART.md) - so this whole file is
skipped unless the image was built with --build-arg ART_SOURCE=....
"""

import pytest

runART = pytest.importorskip(
    "runART",
    reason="ART is licensed separately and is not installed by default - see build/ART.md",
)


def testGetRequiredConfigRaisesOnMissingKeys():
    with pytest.raises(ValueError, match="jobName"):
        runART.getRequiredConfig({"other": 1}, ["jobName", "other"])


def testGetRequiredConfigPassesWhenAllKeysPresent():
    runART.getRequiredConfig({"jobName": "demo", "other": 1}, ["jobName", "other"])


def testBuildPathsDefaultsInputTrainingCsv(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config = {
        "baseDataDir": str(tmp_path),
        "outputDir": "out",
        "jobName": "demo",
        "trainingCsvSuffix": "_train.csv",
        "artReadyCsvSuffix": "_ready.csv",
        "artModelFileName": "model.pkl",
    }
    paths = runART.buildPaths(config)
    assert paths["inputTrainingCsv"] == str(tmp_path / "demo_train.csv")
    assert paths["artReadyCsv"].endswith("demo_ready.csv")
    assert paths["artModelFile"].endswith("model.pkl")


def testBuildPathsHonorsExplicitOverrides(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config = {
        "baseDataDir": str(tmp_path),
        "outputDir": "out",
        "jobName": "demo",
        "trainingCsvSuffix": "_train.csv",
        "artReadyCsvSuffix": "_ready.csv",
        "artModelFileName": "model.pkl",
        "inputTrainingCsv": "custom_input.csv",
    }
    paths = runART.buildPaths(config)
    assert paths["inputTrainingCsv"] == "custom_input.csv"
