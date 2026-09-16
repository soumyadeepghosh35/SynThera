#!/usr/bin/env python
"""Scaffold and run one DORAnet job per starter molecule.

For every SMILES in the input CSV this script creates an individual starter
directory under the output directory, writes a self-contained ``runDORAnet.py``
and ``doranet_config.yaml`` into it, and then runs

    python runDORAnet.py doranet_config.yaml

inside that directory. The per-starter jobs run in parallel across a pool of
worker processes, so each DORAnet run stays isolated in its own directory with
its own configuration and outputs.

Usage:
    python runDORAnet_parallel.py doranet_config_parallel.yaml
"""

import sys
import os
import argparse
import time
import csv
import subprocess
from pathlib import Path
from multiprocessing import Process, Queue, cpu_count

import yaml


# Helper / small molecules excluded from the target set, written into every
# per-starter config so the single-starter runner does not need its own copy.
HELPERS = [
    "O", "O=O", "[H][H]", "O=C=O", "C=O", "[C-]#[O+]", "Br", "[Br][Br]",
    "CO", "C=C", "O=S(O)O", "N", "O=S(=O)(O)O", "O=NO", "N#N",
    "O=[N+]([O-])O", "NO", "C#N", "S", "O=S=O", "N#CO", "[H+]", "OO",
    "Cl", "I", "O=C(O)O", "O=P(O)(O)O", "O=P(O)(O)OP(=O)(O)O", "C",
    "CC", "CC=O", "CC(=O)O", "CCC(=O)O",
]


# The single-starter runner written verbatim into each starter directory. It reads
# only the doranet_config.yaml that sits beside it, so every directory is a complete,
# independently runnable DORAnet job.
RUN_DORANET_TEMPLATE = '''#!/usr/bin/env python
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
'''


def loadConfig(configPath):
    with open(configPath, "r") as handle:
        config = yaml.safe_load(handle)

    defaults = {
        "input": {
            "SMILESfile": None,
            "smiles_column": "Canonical_SMILES",
            "start_index": 0,
            "num_smiles": None,
        },
        "output": {
            "directory": "doranet_output",
        },
        "parallel": {
            "num_workers": 4,
        },
        "network": {
            "generations": 2,
            "total_generations": 2,
            "ruleset": "JN3604IMT",
            "direction": "forward",
        },
        "max_atoms": {
            "C": 12,
            "N": 3,
            "O": 5,
            "S": 3,
        },
        "doranet": {
            "path": None,
        },
        "execution": {
            "run": True,
            "python": sys.executable,
            "timeout_seconds": 3600,
        },
    }

    for section, values in defaults.items():
        if section not in config or config[section] is None:
            config[section] = values
        elif isinstance(values, dict):
            for key, defaultVal in values.items():
                if key not in config[section]:
                    config[section][key] = defaultVal

    return config


def validateConfig(config):
    errors = []

    if not config["input"]["SMILESfile"]:
        errors.append("input.SMILESfile is required")
    elif not os.path.exists(config["input"]["SMILESfile"]):
        errors.append(f"Input file not found: {config['input']['SMILESfile']}")

    if not config["doranet"]["path"]:
        errors.append("doranet.path is required")
    elif not os.path.exists(config["doranet"]["path"]):
        errors.append(f"DORAnet checkout not found: {config['doranet']['path']}")

    if config["parallel"]["num_workers"] < 1:
        errors.append("parallel.num_workers must be at least 1")

    if config["network"]["generations"] < 1:
        errors.append("network.generations must be at least 1")

    if errors:
        print("Configuration errors:")
        for err in errors:
            print(f"  - {err}")
        sys.exit(1)

    config["parallel"]["num_workers"] = min(
        config["parallel"]["num_workers"],
        cpu_count(),
    )

    return config


def printConfig(config):
    print("Configuration")
    print(f"  Input file       : {config['input']['SMILESfile']}")
    print(f"  SMILES column    : {config['input']['smiles_column']}")
    print(f"  Start index      : {config['input']['start_index']}")
    print(f"  Number to process: {config['input']['num_smiles'] or 'all'}")
    print(f"  Output directory : {config['output']['directory']}")
    print(f"  Workers          : {config['parallel']['num_workers']}")
    print(f"  Generations      : {config['network']['generations']}")
    print(f"  Total generations: {config['network']['total_generations']}")
    print(f"  Ruleset          : {config['network']['ruleset']}")
    print(f"  Direction        : {config['network']['direction']}")
    print(f"  Max atoms        : {config['max_atoms']}")
    print(f"  DORAnet path     : {config['doranet']['path']}")
    print(f"  Run per-starter  : {config['execution']['run']}")


def readSmilesFromCsv(csvPath, smilesColumn, startIdx, numSmiles):
    smilesList = []
    with open(csvPath, "r") as handle:
        reader = csv.DictReader(handle)
        if smilesColumn not in reader.fieldnames:
            available = ", ".join(reader.fieldnames)
            raise ValueError(f"Column '{smilesColumn}' not found. Available: {available}")
        for idx, row in enumerate(reader):
            if idx < startIdx:
                continue
            if numSmiles is not None and len(smilesList) >= numSmiles:
                break
            smi = row[smilesColumn]
            if smi and smi.strip():
                smilesList.append(smi.strip())
    return smilesList


def buildStarterConfig(starterSmiles, jobName, config):
    """Assemble the single-starter doranet_config.yaml contents for one starter."""
    return {
        "jobName": jobName,
        "starters": [starterSmiles],
        "doranetPath": config["doranet"]["path"],
        "generations": config["network"]["generations"],
        "totalGenerations": config["network"]["total_generations"],
        "direction": config["network"]["direction"],
        "ruleset": config["network"]["ruleset"],
        "maxAtoms": config["max_atoms"],
        "helpers": config.get("helpers", HELPERS),
    }


def writeStarterFiles(starterIdx, starterSmiles, config):
    """Create the starter directory and write its runDORAnet.py and doranet_config.yaml.
    Returns the directory path and job name."""
    jobName = f"starter_{starterIdx:05d}"
    starterDir = Path(config["output"]["directory"]) / jobName
    starterDir.mkdir(parents=True, exist_ok=True)

    runnerPath = starterDir / "runDORAnet.py"
    runnerPath.write_text(RUN_DORANET_TEMPLATE)

    starterConfig = buildStarterConfig(starterSmiles, jobName, config)
    configPath = starterDir / "doranet_config.yaml"
    with open(configPath, "w") as handle:
        yaml.dump(starterConfig, handle, default_flow_style=False, sort_keys=False)

    return starterDir, jobName


def runOneStarter(starterIdx, starterSmiles, config):
    """Scaffold one starter directory and, if enabled, run its DORAnet job."""
    result = {
        "starter_idx": starterIdx,
        "starter_smiles": starterSmiles,
        "status": "failed",
        "time_seconds": 0,
        "error": None,
    }
    startTime = time.time()

    try:
        starterDir, jobName = writeStarterFiles(starterIdx, starterSmiles, config)

        if not config["execution"]["run"]:
            result["status"] = "scaffolded"
        else:
            completed = subprocess.run(
                [config["execution"]["python"], "runDORAnet.py", "doranet_config.yaml"],
                cwd=str(starterDir),
                capture_output=True,
                text=True,
                timeout=config["execution"]["timeout_seconds"],
            )
            if completed.returncode == 0:
                result["status"] = "success"
            else:
                result["status"] = "failed"
                errorTail = (completed.stderr or "").strip().splitlines()
                result["error"] = errorTail[-1] if errorTail else f"exit code {completed.returncode}"

    except subprocess.TimeoutExpired:
        result["status"] = "failed"
        result["error"] = "timeout"
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = str(exc)

    result["time_seconds"] = round(time.time() - startTime, 2)
    return result


def workerProcess(taskQueue, resultQueue, config, totalTasks):
    """Continuously pull tasks from the queue until the poison pill is received."""
    while True:
        task = taskQueue.get()
        if task is None:
            break

        starterIdx, starterSmiles = task
        result = runOneStarter(starterIdx, starterSmiles, config)

        print(f"[{starterIdx + 1:05d}/{totalTasks:05d}] {result['status'].upper()} | "
              f"Time: {result['time_seconds']}s | "
              f"SMILES: {starterSmiles[:40]}"
              + (f" | Error: {result['error']}" if result["error"] else ""))

        resultQueue.put(result)


def runParallelProcessing(smilesList, config):
    numWorkers = config["parallel"]["num_workers"]
    totalTasks = len(smilesList)

    print(f"\nStarting parallel processing of {totalTasks} starters on {numWorkers} workers")

    taskQueue = Queue()
    for idx, smi in enumerate(smilesList):
        taskQueue.put((idx, smi))
    for _ in range(numWorkers):
        taskQueue.put(None)

    resultQueue = Queue()

    workers = []
    for _ in range(numWorkers):
        p = Process(target=workerProcess, args=(taskQueue, resultQueue, config, totalTasks))
        p.daemon = False
        workers.append(p)
        p.start()

    perResultTimeout = config["execution"]["timeout_seconds"] + 120
    results = []
    for _ in range(totalTasks):
        try:
            results.append(resultQueue.get(timeout=perResultTimeout))
        except Exception as exc:
            print(f"Warning: timeout or error collecting a result: {exc}")
            break

    for p in workers:
        p.join(timeout=60)
        if p.is_alive():
            print(f"Warning: worker {p.pid} did not exit cleanly, terminating")
            p.terminate()
            p.join(timeout=10)

    results.sort(key=lambda r: r["starter_idx"])
    print(f"\nCollected {len(results)} / {totalTasks} results")
    return results


def scaffoldOnly(smilesList, config):
    """Write every starter directory and its files without running DORAnet."""
    results = []
    for idx, smi in enumerate(smilesList):
        result = runOneStarter(idx, smi, config)
        print(f"[{idx + 1:05d}/{len(smilesList):05d}] {result['status'].upper()} | "
              f"SMILES: {smi[:40]}")
        results.append(result)
    return results


def saveSummary(results, outputDir):
    summaryPath = Path(outputDir) / "processing_summary.csv"
    with open(summaryPath, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["starter_idx", "starter_smiles", "status", "time_seconds", "error"])
        for r in results:
            writer.writerow([
                r["starter_idx"],
                r["starter_smiles"],
                r["status"],
                r["time_seconds"],
                r["error"] or "",
            ])
    return summaryPath


def printSummary(results, totalTime):
    successful = [r for r in results if r["status"] in ("success", "scaffolded")]
    failed = [r for r in results if r["status"] == "failed"]
    avgTime = sum(r["time_seconds"] for r in results) / len(results) if results else 0

    print("\nProcessing summary")
    print(f"  Starters processed      : {len(results)}")
    print(f"  Successful / scaffolded : {len(successful)}")
    print(f"  Failed                  : {len(failed)}")
    print(f"  Average time per starter: {avgTime:.2f}s")
    print(f"  Total wall time         : {totalTime:.2f}s")

    if failed:
        print("\nFailed starters:")
        for r in failed[:10]:
            print(f"  [{r['starter_idx']}] {r['starter_smiles'][:50]} | Error: {r['error']}")
        if len(failed) > 10:
            print(f"  ... and {len(failed) - 10} more")


def saveConfigCopy(config, outputDir):
    configCopyPath = Path(outputDir) / "config_used.yaml"
    with open(configCopyPath, "w") as handle:
        yaml.dump(config, handle, default_flow_style=False, sort_keys=False)
    return configCopyPath


def main():
    parser = argparse.ArgumentParser(
        description="Scaffold and run one DORAnet job per starter using a YAML configuration"
    )
    parser.add_argument(
        "config", nargs="?", default="doranet_config_parallel.yaml",
        help="Path to the parallel configuration YAML file",
    )
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: configuration file '{args.config}' not found")
        sys.exit(1)

    print(f"Loading configuration from: {args.config}")
    config = loadConfig(args.config)
    config = validateConfig(config)
    printConfig(config)

    outputDir = Path(config["output"]["directory"])
    outputDir.mkdir(parents=True, exist_ok=True)

    configCopy = saveConfigCopy(config, outputDir)
    print(f"Configuration saved to: {configCopy}")

    print("\nReading SMILES from CSV")
    smilesList = readSmilesFromCsv(
        config["input"]["SMILESfile"],
        config["input"]["smiles_column"],
        config["input"]["start_index"],
        config["input"]["num_smiles"],
    )
    if not smilesList:
        print("Error: no valid SMILES found in input file")
        sys.exit(1)
    print(f"Loaded {len(smilesList)} SMILES")

    totalStartTime = time.time()
    if config["execution"]["run"]:
        results = runParallelProcessing(smilesList, config)
    else:
        print("\nexecution.run is false: scaffolding directories only, not running DORAnet")
        results = scaffoldOnly(smilesList, config)
    totalTime = time.time() - totalStartTime

    summaryPath = saveSummary(results, outputDir)
    print(f"\nSummary saved to: {summaryPath}")

    printSummary(results, totalTime)


if __name__ == "__main__":
    main()
