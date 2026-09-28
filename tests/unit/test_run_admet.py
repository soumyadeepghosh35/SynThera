"""Unit tests for the pure-logic pieces of executables/runADMET.py.

admet_ai is only imported lazily inside loadAdmetModel(), so this module -
and these tests - import fine in either environment without admet_ai
installed.
"""

import pandas as pd
import pytest

import runADMET


def testValidateConfigRaisesOnMissingKeys():
    with pytest.raises(KeyError, match="inputFile"):
        runADMET.validateConfig({"SMILES_Column": "SMILES", "batchSize": 10})


def testValidateConfigAcceptsCompleteConfig():
    runADMET.validateConfig(
        {"inputFile": "x.csv", "SMILES_Column": "SMILES", "batchSize": 10}
    )


def testBuildPathsRaisesWhenInputMissing(tmp_path):
    with pytest.raises(FileNotFoundError):
        runADMET.buildPaths({"inputFile": str(tmp_path / "missing.csv")})


def testBuildPathsUsesDefaultSuffixes(tmp_path):
    inputCsv = tmp_path / "molecules.csv"
    inputCsv.write_text("Canonical_SMILES\nCCO\n")

    paths = runADMET.buildPaths({"inputFile": str(inputCsv), "outputDir": str(tmp_path)})

    assert paths["outputCsvPath"] == str(tmp_path / "molecules_wADMET.csv")
    assert paths["diversityCsvPath"] == str(tmp_path / "molecules_diverseInput.csv")
    assert paths["inputFile"] == "molecules.csv"


def testBuildPathsHonorsCustomSuffixes(tmp_path):
    inputCsv = tmp_path / "molecules.csv"
    inputCsv.write_text("Canonical_SMILES\nCCO\n")

    paths = runADMET.buildPaths(
        {
            "inputFile": str(inputCsv),
            "outputDir": str(tmp_path),
            "outputSuffix": "_scored",
            "diversityOutputSuffix": "_diverse",
        }
    )

    assert paths["outputCsvPath"] == str(tmp_path / "molecules_scored.csv")
    assert paths["diversityCsvPath"] == str(tmp_path / "molecules_diverse.csv")


def testCanonicalizeSmilesReturnsNoneForInvalidInput():
    assert runADMET.canonicalizeSmiles("not_a_real_smiles") is None


def testCanonicalizeSmilesReturnsNoneForNan():
    assert runADMET.canonicalizeSmiles(float("nan")) is None


def testCanonicalizeSmilesCanonicalizesEquivalentInput():
    # Two spellings of ethanol should canonicalize to the same SMILES.
    a = runADMET.canonicalizeSmiles("CCO")
    b = runADMET.canonicalizeSmiles("OCC")
    assert a == b


def testDeduplicateForPredictionDropsDuplicatesAndNulls():
    df = pd.DataFrame({"SMILES": ["CCO", "CCO", None, "CCN"]})
    deduped = runADMET.deduplicateForPrediction(df, "SMILES")
    assert sorted(deduped["SMILES"]) == ["CCN", "CCO"]
