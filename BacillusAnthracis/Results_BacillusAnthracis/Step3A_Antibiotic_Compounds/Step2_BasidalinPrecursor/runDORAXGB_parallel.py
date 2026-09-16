#!/usr/bin/env python
"""Score every starter's reaction network with DORA-XGB, in process, one pool of workers.

For each starter_* directory produced by the DORAnet run, this reads (or builds from
that directory's *_network_pretreated.json) a reaction table and writes
reactionDF_wDORAXGBfeasibility.csv into the directory, which is the file the Step3
aggregation reads.

Unlike a per-directory subprocess design, this starts a small fixed pool of worker
processes, loads the four DORA-XGB models once per worker, and reuses them for every
directory that worker handles. That removes thousands of interpreter startups and
model reloads, which is what keeps the machine's load average down.

Usage:
    python runDORAXGB_parallel.py doraxgb_config_parallel.yaml
"""

import os

# Thread limits MUST be set before DORA-XGB / numpy import, in the parent and in every
# spawned worker (each worker re-imports this module, so these run there too).
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["BLIS_NUM_THREADS"] = "1"
os.environ["OMP_DYNAMIC"] = "FALSE"
os.environ["MKL_DYNAMIC"] = "FALSE"

import sys
import argparse
import time
import csv
import json
import glob
import multiprocessing as mp
from pathlib import Path
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor, as_completed

import yaml
import pandas as pd
from tqdm import tqdm

from DORA_XGB import DORA_XGB


# All four DORA-XGB cofactor-positioning rules, always scored together.
RULE_TO_POSITIONING = {
    "rule1": "by_descending_MW",
    "rule2": "by_ascending_MW",
    "rule3": "add_concat",
    "rule4": "add_subtract",
}
RULES = list(RULE_TO_POSITIONING.keys())

# Per-worker state, populated once by initWorker and reused for every directory.
_WORKER = {}


def initWorker(params):
    """Runs once per worker process: load the four models a single time and build a
    per-worker reaction cache. Every directory this worker handles reuses them."""
    models = {
        r: DORA_XGB.feasibility_classifier(cofactor_positioning=RULE_TO_POSITIONING[r])
        for r in RULES
    }

    @lru_cache(maxsize=int(params.get("cacheSize", 200000)))
    def predictReaction(rxnStr):
        scores = {}
        for r in RULES:
            model = models[r]
            scores[f"feasibilityScore_{r}"] = model.predict_proba(rxnStr)
            scores[f"feasibilityLabel_{r}"] = model.predict_label(rxnStr)
        return scores

    _WORKER["predict"] = predictReaction
    _WORKER["params"] = params


def buildReactionDf(jsonPath):
    """Parse a DORAnet *_network_pretreated.json into a reaction DataFrame with a
    reactionString column."""
    with open(jsonPath, "r", encoding="utf-8") as handle:
        reactionList = json.load(handle)

    rows = []
    for rxn in reactionList:
        parts = str(rxn).split(">")
        if len(parts) != 4:
            continue
        reactants, ruleName, _metaBlock, products = parts
        rows.append({
            "reactants": reactants.strip(),
            "products": products.strip(),
            "reactionString": f"{reactants.strip()} >> {products.strip()}",
            "ruleName": ruleName.strip(),
        })
    return pd.DataFrame(rows)


def loadOrBuildReactionDf(dirPath, params):
    """Read the directory's reaction CSV, or build it from the network JSON and cache it."""
    reactionCsvPath = dirPath / params["reactionCSV"]
    if reactionCsvPath.exists():
        return pd.read_csv(reactionCsvPath)

    jsonMatches = sorted(glob.glob(str(dirPath / params["networkJsonGlob"])))
    if not jsonMatches:
        raise FileNotFoundError(
            f"no {params['reactionCSV']} and no network json to build it from"
        )
    df = buildReactionDf(jsonMatches[0])
    df.to_csv(reactionCsvPath, index=False)
    return df


def scoreDirectory(dirPath):
    """Score one starter directory in process, reusing this worker's loaded models."""
    result = {"dir": os.path.basename(dirPath), "status": "failed",
              "n_reactions": -1, "time_seconds": 0, "error": None}
    startTime = time.time()

    try:
        params = _WORKER["params"]
        predict = _WORKER["predict"]
        dirPath = Path(dirPath)

        df = loadOrBuildReactionDf(dirPath, params)
        if "reactionString" not in df.columns:
            raise ValueError("'reactionString' column not found")
        if params.get("testMode", False):
            df = df.iloc[:int(params.get("testRows", 100))].copy()

        df["rxn_str"] = df["reactionString"].astype(str).str.replace(" ", "", regex=False)
        uniqueRxns = pd.Series(df["rxn_str"].dropna().unique()).tolist()

        if uniqueRxns:
            scored = pd.DataFrame([{**predict(rxn), "rxn_str": rxn} for rxn in uniqueRxns])
            out = df.merge(scored, on="rxn_str", how="left")
        else:
            out = df

        out.to_csv(dirPath / "reactionDF_wDORAXGBfeasibility.csv", index=False)
        result["n_reactions"] = len(df)
        result["status"] = "success"

    except Exception as exc:
        result["error"] = str(exc)

    result["time_seconds"] = round(time.time() - startTime, 2)
    return result


def loadConfig(configPath):
    with open(configPath, "r") as handle:
        config = yaml.safe_load(handle)

    defaults = {
        "input": {
            "outputDir": "doranet_output",
            "starterGlob": "starter_*",
            "reactionCSV": "reactionDF.csv",
            "networkJsonGlob": "*_network_pretreated.json",
        },
        "parallel": {
            "num_workers": 4,
        },
        "scoring": {
            "testMode": False,
            "testRows": 100,
            "cacheSize": 200000,
        },
        "execution": {
            "run": True,
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
    if not os.path.isdir(config["input"]["outputDir"]):
        errors.append(f"input.outputDir not found: {config['input']['outputDir']}")
    if config["parallel"]["num_workers"] < 1:
        errors.append("parallel.num_workers must be at least 1")

    if errors:
        print("Configuration errors:")
        for err in errors:
            print(f"  - {err}")
        sys.exit(1)

    config["parallel"]["num_workers"] = min(config["parallel"]["num_workers"], mp.cpu_count())
    return config


def printConfig(config, nDirs):
    print("Configuration")
    print(f"  Output directory  : {config['input']['outputDir']}")
    print(f"  Starter glob      : {config['input']['starterGlob']}")
    print(f"  Starter dirs found: {nDirs}")
    print(f"  Reaction CSV      : {config['input']['reactionCSV']}")
    print(f"  Workers           : {config['parallel']['num_workers']}  (models loaded once per worker)")
    print(f"  Cache size        : {config['scoring']['cacheSize']}")
    print(f"  Run scoring       : {config['execution']['run']}")


def discoverStarterDirs(config):
    root = config["input"]["outputDir"]
    return sorted(
        p for p in glob.glob(os.path.join(root, config["input"]["starterGlob"]))
        if os.path.isdir(p)
    )


def workerParams(config):
    return {
        "reactionCSV": config["input"]["reactionCSV"],
        "networkJsonGlob": config["input"]["networkJsonGlob"],
        "testMode": config["scoring"]["testMode"],
        "testRows": config["scoring"]["testRows"],
        "cacheSize": config["scoring"]["cacheSize"],
    }


def runInProcess(starterDirs, config):
    numWorkers = config["parallel"]["num_workers"]
    params = workerParams(config)
    print(f"\nScoring {len(starterDirs)} directories on {numWorkers} workers "
          f"(each worker loads the models once)")

    results = []
    nSuccess = nFailed = 0
    # spawn keeps each worker's model load clean and independent; the models are loaded
    # inside initWorker, so nothing heavy is pickled across the process boundary.
    spawnContext = mp.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=numWorkers,
        mp_context=spawnContext,
        initializer=initWorker,
        initargs=(params,),
    ) as executor:
        futures = {executor.submit(scoreDirectory, d): d for d in starterDirs}

        progressBar = tqdm(as_completed(futures), total=len(futures),
                           desc="Scoring directories", unit="dir")
        for future in progressBar:
            result = future.result()
            results.append(result)
            if result["status"] == "success":
                nSuccess += 1
            else:
                nFailed += 1
                progressBar.write(f"  FAILED {result['dir']}: {result['error']}")
            progressBar.set_postfix(ok=nSuccess, failed=nFailed, last=result["dir"], refresh=False)

    results.sort(key=lambda r: r["dir"])
    return results


def scaffoldOnly(starterDirs, config):
    """Build each directory's reaction CSV from its network JSON without scoring."""
    params = workerParams(config)
    results = []
    for dirPath in tqdm(starterDirs, desc="Building reaction CSVs", unit="dir"):
        result = {"dir": os.path.basename(dirPath), "status": "failed",
                  "n_reactions": -1, "time_seconds": 0, "error": None}
        startTime = time.time()
        try:
            df = loadOrBuildReactionDf(Path(dirPath), params)
            result["n_reactions"] = len(df)
            result["status"] = "scaffolded"
        except Exception as exc:
            result["error"] = str(exc)
        result["time_seconds"] = round(time.time() - startTime, 2)
        results.append(result)
    return results


def saveSummary(results, outputDir):
    summaryPath = Path(outputDir) / "doraxgb_processing_summary.csv"
    with open(summaryPath, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["dir", "status", "n_reactions", "time_seconds", "error"])
        for r in results:
            writer.writerow([r["dir"], r["status"], r["n_reactions"], r["time_seconds"], r["error"] or ""])
    return summaryPath


def printSummary(results, totalTime):
    successful = [r for r in results if r["status"] in ("success", "scaffolded")]
    failed = [r for r in results if r["status"] == "failed"]
    avgTime = sum(r["time_seconds"] for r in results) / len(results) if results else 0

    print("\nProcessing summary")
    print(f"  Directories processed   : {len(results)}")
    print(f"  Successful / scaffolded : {len(successful)}")
    print(f"  Failed                  : {len(failed)}")
    print(f"  Average time per dir    : {avgTime:.2f}s")
    print(f"  Total wall time         : {totalTime:.2f}s")

    if failed:
        print("\nFailed directories:")
        for r in failed[:10]:
            print(f"  {r['dir']} | Error: {r['error']}")
        if len(failed) > 10:
            print(f"  ... and {len(failed) - 10} more")


def saveConfigCopy(config, outputDir):
    configCopyPath = Path(outputDir) / "doraxgb_config_used.yaml"
    with open(configCopyPath, "w") as handle:
        yaml.dump(config, handle, default_flow_style=False, sort_keys=False)
    return configCopyPath


def main():
    parser = argparse.ArgumentParser(
        description="Score every starter directory with DORA-XGB using an in-process worker pool"
    )
    parser.add_argument("config", nargs="?", default="doraxgb_config_parallel.yaml",
                        help="Path to the parallel configuration YAML file")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: configuration file '{args.config}' not found")
        sys.exit(1)

    print(f"Loading configuration from: {args.config}")
    config = loadConfig(args.config)
    config = validateConfig(config)

    starterDirs = discoverStarterDirs(config)
    printConfig(config, len(starterDirs))
    if not starterDirs:
        print("Error: no starter directories matched; nothing to do")
        sys.exit(1)

    configCopy = saveConfigCopy(config, config["input"]["outputDir"])
    print(f"Configuration saved to: {configCopy}")

    totalStartTime = time.time()
    if config["execution"]["run"]:
        results = runInProcess(starterDirs, config)
    else:
        print("\nexecution.run is false: building reaction CSVs only, not scoring")
        results = scaffoldOnly(starterDirs, config)
    totalTime = time.time() - totalStartTime

    summaryPath = saveSummary(results, config["input"]["outputDir"])
    print(f"\nSummary saved to: {summaryPath}")
    printSummary(results, totalTime)


if __name__ == "__main__":
    main()
