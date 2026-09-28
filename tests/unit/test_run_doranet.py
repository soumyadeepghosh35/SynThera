"""Unit tests for the config-handling logic in executables/runDORAnet.py.

Only loadConfig/resolveConfig are exercised here - pure functions that need
nothing beyond PyYAML, so this file runs the same way in either environment.
"""

import pytest

import runDORAnet


def writeYaml(tmpPath, text):
    configPath = tmpPath / "config.yaml"
    configPath.write_text(text)
    return configPath


def testLoadConfigMissingFileRaises(tmp_path):
    with pytest.raises(FileNotFoundError):
        runDORAnet.loadConfig(tmp_path / "does_not_exist.yaml")


def testLoadConfigReadsYaml(tmp_path):
    configPath = writeYaml(
        tmp_path, "jobName: demo\nstarters:\n  - CCO\ndoranetPath: /tmp/doranet\n"
    )
    config = runDORAnet.loadConfig(configPath)
    assert config == {
        "jobName": "demo",
        "starters": ["CCO"],
        "doranetPath": "/tmp/doranet",
    }


def testResolveConfigMissingRequiredKeysRaises():
    with pytest.raises(ValueError, match="jobName"):
        runDORAnet.resolveConfig({"starters": ["CCO"], "doranetPath": "/tmp/doranet"})


def testResolveConfigFillsDefaults():
    settings = runDORAnet.resolveConfig(
        {"jobName": "demo", "starters": ["CCO"], "doranetPath": "/tmp/doranet"}
    )
    assert settings["generations"] == 3
    assert settings["totalGenerations"] == 3
    assert settings["direction"] == "forward"
    assert settings["ruleset"] == "JN3604IMT"
    assert settings["maxAtoms"] == runDORAnet.DEFAULT_MAX_ATOMS
    assert settings["helpers"] == set(runDORAnet.DEFAULT_HELPERS)


def testResolveConfigCoercesSingleStarterStringToSet():
    settings = runDORAnet.resolveConfig(
        {"jobName": "demo", "starters": "CCO", "doranetPath": "/tmp/doranet"}
    )
    assert settings["starters"] == {"CCO"}


def testResolveConfigHonorsExplicitOverrides():
    settings = runDORAnet.resolveConfig(
        {
            "jobName": "demo",
            "starters": ["CCO", "CCN"],
            "doranetPath": "/tmp/doranet",
            "helpers": ["O"],
            "generations": 2,
            "totalGenerations": 5,
            "direction": "retro",
            "ruleset": "custom_ruleset",
            "maxAtoms": {"C": 4},
        }
    )
    assert settings["helpers"] == {"O"}
    assert settings["generations"] == 2
    assert settings["totalGenerations"] == 5
    assert settings["direction"] == "retro"
    assert settings["ruleset"] == "custom_ruleset"
    assert settings["maxAtoms"] == {"C": 4}
