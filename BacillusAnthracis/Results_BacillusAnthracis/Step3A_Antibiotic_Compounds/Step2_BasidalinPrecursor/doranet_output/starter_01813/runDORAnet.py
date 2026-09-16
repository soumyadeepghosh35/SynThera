#!/usr/bin/env python
"""Single-starter DORAnet runner.

Reads the doranet_config.yaml in the current directory and generates the reaction
network, the molecules table, and the pathway files for one starter molecule.

Invoked by runDORAnet_parallel.py, or directly as:
    python runDORAnet.py doranet_config.yaml
"""

import sys
import argparse
import csv

import yaml
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def loadConfig(configPath):
    with open(configPath, "r") as handle:
        return yaml.safe_load(handle)


def writeMoleculesCsv(network, starters, jobName):
    smilesList = [mol.uid for mol in network.mols]
    outputPath = f"{jobName}_molecules.csv"
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
    return smilesList


def main():
    parser = argparse.ArgumentParser(description="Single-starter DORAnet runner")
    parser.add_argument(
        "config", nargs="?", default="doranet_config.yaml",
        help="Path to this starter's doranet_config.yaml",
    )
    args = parser.parse_args()

    config = loadConfig(args.config)

    doranetPath = config["doranetPath"]
    sys.path.insert(0, str(doranetPath))
    import doranet.modules.enzymatic as enzymatic
    import doranet.modules.post_processing as post_processing

    jobName = config["jobName"]
    starters = set(config["starters"])
    helpers = set(config.get("helpers", []))
    generations = int(config["generations"])
    totalGenerations = int(config.get("totalGenerations", generations))
    maxAtoms = config["maxAtoms"]
    direction = config.get("direction", "forward")
    ruleset = config.get("ruleset", "JN3604IMT")

    forwardNetwork = enzymatic.generate_network(
        job_name=jobName,
        starters=starters,
        gen=generations,
        max_atoms=maxAtoms,
        direction=direction,
        ruleset=ruleset,
    )

    smilesList = writeMoleculesCsv(forwardNetwork, starters, jobName)
    allTargets = set(smilesList) - starters - helpers

    if allTargets:
        post_processing.one_step(
            networks={forwardNetwork},
            total_generations=totalGenerations,
            starters=starters,
            helpers=helpers,
            target=allTargets,
            job_name=jobName,
        )

    print(f"{jobName}: {len(smilesList)} molecules, {len(allTargets)} targets")


if __name__ == "__main__":
    main()
