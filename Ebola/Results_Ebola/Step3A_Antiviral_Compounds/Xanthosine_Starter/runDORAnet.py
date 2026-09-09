#!/usr/bin/env python
"""Run DORAnet enzymatic network generation and pathway search from a config file.

Usage:
    python runDORAnet.py doranet_config.yaml

Only the YAML file changes between runs, so this script is reusable across
different starters and rulesets.
"""

import sys
import csv
import time
from pathlib import Path

import yaml

# Default cofactor/helper molecules, used when the config omits 'helpers'.
DEFAULT_HELPERS = [
    "O", "O=O", "[H][H]", "O=C=O", "C=O", "[C-]#[O+]", "Br", "[Br][Br]", "CO",
    "C=C", "O=S(O)O", "N", "O=S(=O)(O)O", "O=NO", "N#N", "O=[N+]([O-])O", "NO",
    "C#N", "S", "O=S=O", "N#CO",
]

DEFAULT_MAX_ATOMS = {"C": 15, "N": 6, "O": 8, "S": 3}


def loadConfig(configPath):
    """Read the YAML config file into a dictionary."""
    path = Path(configPath)
    if not path.is_file():
        raise FileNotFoundError(f"Config file not found: {path}")
    with open(path) as handle:
        config = yaml.safe_load(handle) or {}
    return config


def resolveConfig(config):
    """Validate required fields, fill defaults, and return a clean settings dict."""
    required = ["jobName", "starters", "doranetPath"]
    missing = [key for key in required if not config.get(key)]
    if missing:
        raise ValueError(f"Missing required config keys: {missing}")

    starters = config["starters"]
    if isinstance(starters, str):
        starters = [starters]

    helpers = config.get("helpers") or DEFAULT_HELPERS
    generations = int(config.get("generations", 3))

    return {
        "doranetPath": str(config["doranetPath"]),
        "jobName": str(config["jobName"]),
        "starters": set(starters),
        "helpers": set(helpers),
        "generations": generations,
        "direction": str(config.get("direction", "forward")),
        "ruleset": str(config.get("ruleset", "JN3604IMT")),
        "maxAtoms": dict(config.get("maxAtoms", DEFAULT_MAX_ATOMS)),
        "totalGenerations": int(config.get("totalGenerations", generations)),
    }


def writeMoleculeCsv(smilesList, starters, outputPath):
    """Write one row per molecule with formula, exact mass, and heavy-atom count."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    with open(outputPath, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["SMILES", "Is_Starter", "MolFormula", "MolWeight", "NumHeavyAtoms"])
        for smi in smilesList:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                formula = rdMolDescriptors.CalcMolFormula(mol)
                molWeight = round(rdMolDescriptors.CalcExactMolWt(mol), 4)
                numHeavy = mol.GetNumHeavyAtoms()
            else:
                formula, molWeight, numHeavy = "N/A", 0, 0
            writer.writerow([smi, smi in starters, formula, molWeight, numHeavy])


def main():
    if len(sys.argv) != 2:
        print("Usage: python runDORAnet.py doranet_config.yaml")
        sys.exit(1)

    startTime = time.time()

    config = loadConfig(sys.argv[1])
    settings = resolveConfig(config)

    doranetPath = Path(settings["doranetPath"])
    if not doranetPath.is_dir():
        raise FileNotFoundError(f"doranetPath does not exist: {doranetPath}")
    sys.path.insert(0, str(doranetPath))

    import doranet.modules.enzymatic as enzymatic
    import doranet.modules.post_processing as post_processing

    starters = settings["starters"]
    helpers = settings["helpers"]
    jobName = settings["jobName"]

    print(f"Job: {jobName}")
    print(f"Starters ({len(starters)}): {sorted(starters)}")

    forwardNetwork = enzymatic.generate_network(
        job_name=jobName,
        starters=starters,
        gen=settings["generations"],
        max_atoms=settings["maxAtoms"],
        direction=settings["direction"],
        ruleset=settings["ruleset"],
    )

    # Include the starters, then every newly generated molecule.
    generatedUids = [mol.uid for mol in forwardNetwork.mols if mol.uid not in starters]
    smilesList = list(starters) + generatedUids
    print(f"Generated {len(generatedUids)} new molecules + {len(starters)} starters")

    # Every generated molecule (excluding starters and helpers) becomes a target.
    allTargets = set(smilesList) - starters - helpers
    print(f"Using {len(allTargets)} molecules as pathway targets")

    outputPath = Path(f"{jobName}_molecules.csv")
    writeMoleculeCsv(smilesList, starters, outputPath)
    print(f"Saved molecules to {outputPath}")

    if allTargets:
        post_processing.one_step(
            networks={forwardNetwork},
            total_generations=settings["totalGenerations"],
            starters=starters,
            helpers=helpers,
            target=allTargets,
            job_name=jobName,
        )
        print(f"Pathway files generated with prefix: {jobName}")
    else:
        print("No targets found for pathway generation")

    print(f"Time: {time.time() - startTime:.2f} s")


if __name__ == "__main__":
    main()
