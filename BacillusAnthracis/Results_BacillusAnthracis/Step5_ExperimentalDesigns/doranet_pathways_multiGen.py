#!/usr/bin/env python3
"""
Run DORAnet pathway jobs from YAML config for multiple generations.

Usage:
    python doranet_pathways_multiGen.py config_doranet_multiGen.yaml

Key behavior:
    1. Reads starter_Canonical_SMILES and Canonical_SMILES (target) from targetsCsvPath.
    2. Computes maxAtoms dynamically per starter: ceil(atomCount * atomMultiplier).
    3. Runs separate jobs for each (starter, target) pair x generation in generationsToRun.
    4. Searches pathways specifically to the user-requested target, not all generated molecules.
    5. Optionally filters DORAnet's <=N-step pathways to exactly N-step pathways.
"""

from __future__ import annotations

import ast
import csv
import math
import os
import shutil
import sys
import time
import textwrap
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def setThreadEnv() -> None:
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"


def loadYaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def resolveGenerationRuns(config: dict[str, Any]) -> list[int]:
    raw = config.get("generationsToRun", config.get("generations", 3))
    if isinstance(raw, int):
        raw = [raw]
    return sorted({int(g) for g in raw if int(g) > 0})


def getSmilesProps(smiles: str) -> tuple[str, float, int]:
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "N/A", 0.0, 0
    return (
        rdMolDescriptors.CalcMolFormula(mol),
        round(rdMolDescriptors.CalcExactMolWt(mol), 4),
        mol.GetNumHeavyAtoms(),
    )


def canonicalSmiles(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(str(smiles))
    return Chem.MolToSmiles(mol, canonical=True) if mol else None


def computeMaxAtoms(starterSmiles: str, multiplier: float = 1.5) -> dict[str, int]:
    """
    Compute maxAtoms for DORAnet expansion from starter atom counts.
    Each atom type count is scaled by multiplier and rounded up (ceil).
    Minimum value is 1 for each type to avoid blocking any chemistry.
    """
    mol = Chem.MolFromSmiles(str(starterSmiles))
    if mol is None:
        raise ValueError(f"Could not parse starter SMILES: {starterSmiles}")
    counts = {"C": 0, "N": 0, "O": 0, "S": 0}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in counts:
            counts[atom.GetSymbol()] += 1
    return {sym: max(1, math.ceil(cnt * multiplier)) for sym, cnt in counts.items()}


def loadStarterTargetPairs(
    csvPath: Path,
    starterCol: str,
    targetCol: str,
    targetsToTest: list[int] | None = None,
) -> list[tuple[str, str, int]]:
    """
    Read (starterSmiles, targetSmiles, originalIndex) from CSV.
    Returns unique pairs in CSV row order, optionally filtered by 1-based row index.
    """
    df = pd.read_csv(csvPath, dtype=str)
    if starterCol not in df.columns:
        raise ValueError(f"Column '{starterCol}' not found in {csvPath}. Available: {list(df.columns)}")
    if targetCol not in df.columns:
        raise ValueError(f"Column '{targetCol}' not found in {csvPath}. Available: {list(df.columns)}")

    df = df[[starterCol, targetCol]].dropna().reset_index(drop=True)
    keepIdx = set(targetsToTest) if targetsToTest else None

    seen: set[tuple[str, str]] = set()
    pairs: list[tuple[str, str, int]] = []
    for rowIdx, row in df.iterrows():
        idx = int(rowIdx) + 1
        if keepIdx is not None and idx not in keepIdx:
            continue
        starter = str(row[starterCol]).strip()
        target  = str(row[targetCol]).strip()
        if not starter or not target:
            continue
        if (starter, target) not in seen:
            seen.add((starter, target))
            pairs.append((starter, target, idx))

    return pairs


def countPathwayBlockSteps(blockLines: list[str]) -> int | None:
    prefix = "reaction SMILES stoichiometry "
    for line in blockLines:
        if line.startswith(prefix):
            payload = line[len(prefix):].strip()
            try:
                parsed = ast.literal_eval(payload)
                if isinstance(parsed, list):
                    return len(parsed)
            except Exception:
                return None
    return None


def splitPathwayBlocks(pathwaysTxtPath: Path) -> list[list[str]]:
    if not pathwaysTxtPath.exists():
        return []
    lines = pathwaysTxtPath.read_text(encoding="utf-8").splitlines()
    blocks: list[list[str]] = []
    current: list[str] = []
    for line in lines:
        if line.startswith("pathway number ") and current:
            blocks.append(current)
            current = [line]
        else:
            current.append(line)
    if current:
        blocks.append(current)
    return blocks


def filterPathwaysTxtByExactSteps(jobName: str, exactSteps: int) -> dict[str, Any]:
    pathwaysPath = Path(f"{jobName}_pathways.txt")
    backupPath   = Path(f"{jobName}_pathways_unfiltered.txt")
    exactPath    = Path(f"{jobName}_pathways_exact{exactSteps}.txt")

    if not pathwaysPath.exists():
        return {
            "pathwaysTxtFound": False,
            "numPathwayBlocksOriginal": 0,
            "numPathwayBlocksExact": 0,
            "exactPathwaysTxt": str(exactPath),
        }

    blocks      = splitPathwayBlocks(pathwaysPath)
    exactBlocks = [b for b in blocks if countPathwayBlockSteps(b) == exactSteps]

    if not backupPath.exists():
        shutil.copy2(pathwaysPath, backupPath)

    exactText = "\n\n".join("\n".join(block) for block in exactBlocks)
    if exactText:
        exactText += "\n"
    exactPath.write_text(exactText, encoding="utf-8")
    pathwaysPath.write_text(exactText, encoding="utf-8")

    return {
        "pathwaysTxtFound": True,
        "numPathwayBlocksOriginal": len(blocks),
        "numPathwayBlocksExact": len(exactBlocks),
        "exactPathwaysTxt": str(exactPath),
        "backupPathwaysTxt": str(backupPath),
    }


def maybeRunPostProcessing(
    network: Any,
    jobName: str,
    userTarget: set[str],
    starters: set[str],
    helpers: set[str],
    config: dict[str, Any],
    generationRun: int,
) -> dict[str, Any]:
    import doranet.modules.post_processing as postProcessing

    searchDepth = generationRun if bool(config.get("searchDepthEqualsGeneration", True)) else int(config.get("searchDepth", generationRun))
    maxNumRxns  = generationRun if bool(config.get("maxNumRxnsEqualsGeneration", True)) else int(config.get("maxNumRxns", generationRun))
    minRxnAtomEconomy   = float(config.get("minRxnAtomEconomy", 0.0))
    runRanking          = bool(config.get("runRanking", True))
    runVisualization    = bool(config.get("runVisualization", False))
    numProcess          = int(config.get("numProcess", 1))

    postProcessing.pretreat_networks(
        networks={network},
        total_generations=generationRun,
        starters=starters,
        helpers=helpers,
        job_name=jobName,
        remove_pure_helpers_rxns=bool(config.get("removePureHelpersRxns", False)),
        sanitize=bool(config.get("sanitize", True)),
        transform_enols_flag=bool(config.get("transformEnolsFlag", False)),
    )

    postProcessing.pathway_finder(
        starters=starters,
        helpers=helpers,
        target=userTarget,
        search_depth=searchDepth,
        max_num_rxns=maxNumRxns,
        min_rxn_atom_economy=minRxnAtomEconomy,
        job_name=jobName,
        consider_name_difference=bool(config.get("considerNameDifference", True)),
    )

    filterInfo: dict[str, Any] = {"pathwaysTxtFound": Path(f"{jobName}_pathways.txt").exists()}

    if bool(config.get("filterExactNumRxns", True)):
        exactSteps = generationRun if bool(config.get("exactNumRxnsEqualsGeneration", True)) else int(config.get("exactNumRxns", generationRun))
        filterInfo = filterPathwaysTxtByExactSteps(jobName, exactSteps)

    if runRanking:
        postProcessing.pathway_ranking(
            starters=starters,
            helpers=helpers,
            target=userTarget,
            weights=config.get("weights", None),
            num_process=numProcess,
            reaxys_result_name=config.get("reaxysResultName", None),
            job_name=jobName,
            cool_reactions=config.get("coolReactions", None),
        )

    if runVisualization:
        postProcessing.pathway_visualization(
            starters=starters,
            helpers=helpers,
            num_process=numProcess,
            reaxys_result_name=config.get("reaxysResultName", "default"),
            job_name=jobName,
            exclude_smiles=config.get("excludeSmiles", None),
            reaxys_rxn_color=config.get("reaxysRxnColor", "blue"),
            normal_rxn_color=config.get("normalRxnColor", "black"),
        )

    return {
        "generationRun": generationRun,
        "searchDepth": searchDepth,
        "maxNumRxns": maxNumRxns,
        "minRxnAtomEconomy": minRxnAtomEconomy,
        "targetUsedForPathwaySearch": ";".join(sorted(userTarget)),
        "runRanking": runRanking,
        "runVisualization": runVisualization,
        **filterInfo,
    }


def saveReproScript(
    jobDir: Path,
    config: dict[str, Any],
    jobName: str,
    starterSmiles: str,
    targetSmiles: str,
    maxAtoms: dict[str, int],
    generationRun: int,
) -> None:
    searchDepth = generationRun if bool(config.get("searchDepthEqualsGeneration", True)) else int(config.get("searchDepth", generationRun))
    maxNumRxns  = generationRun if bool(config.get("maxNumRxnsEqualsGeneration", True)) else int(config.get("maxNumRxns", generationRun))
    exactNumRxns = generationRun if bool(config.get("exactNumRxnsEqualsGeneration", True)) else int(config.get("exactNumRxns", generationRun))

    scriptText = f'''\
#!/usr/bin/env python3
import ast, csv, os, shutil, sys, time
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

doranetPath = Path(r"{config['doranetPath']}")
sys.path.insert(0, str(doranetPath))

import doranet.modules.enzymatic as enzymatic
import doranet.modules.post_processing as postProcessing

jobName      = "{jobName}"
starters     = {repr({starterSmiles})}
helpers      = {repr(set(config['helpers']))}
target       = {repr({targetSmiles})}
maxAtoms     = {repr(maxAtoms)}
generations  = {generationRun}
ruleset      = "{config['ruleset']}"
searchDepth  = {searchDepth}
maxNumRxns   = {maxNumRxns}
minRxnAtomEconomy   = {float(config.get('minRxnAtomEconomy', 0.0))}
filterExactNumRxns  = {bool(config.get('filterExactNumRxns', True))}
exactNumRxns        = {exactNumRxns}
runRanking          = {bool(config.get('runRanking', True))}
runVisualization    = {bool(config.get('runVisualization', False))}

def getSmilesProps(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "N/A", 0.0, 0
    return rdMolDescriptors.CalcMolFormula(mol), round(rdMolDescriptors.CalcExactMolWt(mol), 4), mol.GetNumHeavyAtoms()

def splitPathwayBlocks(p):
    if not p.exists():
        return []
    lines = p.read_text(encoding="utf-8").splitlines()
    blocks, current = [], []
    for line in lines:
        if line.startswith("pathway number ") and current:
            blocks.append(current)
            current = [line]
        else:
            current.append(line)
    if current:
        blocks.append(current)
    return blocks

def countSteps(block):
    import ast
    prefix = "reaction SMILES stoichiometry "
    for line in block:
        if line.startswith(prefix):
            try:
                parsed = ast.literal_eval(line[len(prefix):].strip())
                return len(parsed) if isinstance(parsed, list) else None
            except Exception:
                return None
    return None

t0 = time.time()
network = enzymatic.generate_network(
    job_name=jobName, starters=starters, gen=generations, max_atoms=maxAtoms,
    direction="{config.get('direction', 'forward')}", targets=target, ruleset=ruleset,
)

smilesList = list(starters) + [m.uid for m in network.mols if m.uid not in starters]
with Path(f"{{jobName}}_molecules.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["SMILES", "isStarter", "molFormula", "molWeight", "numHeavyAtoms"])
    for s in smilesList:
        w.writerow([s, s in starters, *getSmilesProps(s)])

postProcessing.pretreat_networks(networks={{network}}, total_generations=generations, starters=starters, helpers=helpers, job_name=jobName)
postProcessing.pathway_finder(starters=starters, helpers=helpers, target=target, search_depth=searchDepth, max_num_rxns=maxNumRxns, min_rxn_atom_economy=minRxnAtomEconomy, job_name=jobName)

if filterExactNumRxns:
    pPath = Path(f"{{jobName}}_pathways.txt")
    ePath  = Path(f"{{jobName}}_pathways_exact{{exactNumRxns}}.txt")
    if pPath.exists():
        blocks = splitPathwayBlocks(pPath)
        exact  = [b for b in blocks if countSteps(b) == exactNumRxns]
        txt    = "\\n\\n".join("\\n".join(b) for b in exact)
        ePath.write_text(txt + "\\n" if txt else "", encoding="utf-8")
        pPath.write_text(txt + "\\n" if txt else "", encoding="utf-8")
        print(f"Exact-step filter Gen{{generations}}: {{len(exact)}}/{{len(blocks)}} pathways retained")

if runRanking:
    postProcessing.pathway_ranking(starters=starters, helpers=helpers, target=target, job_name=jobName, num_process=1)
if runVisualization:
    postProcessing.pathway_visualization(starters=starters, helpers=helpers, job_name=jobName, num_process=1)

print(f"Done {{jobName}} in {{time.time() - t0:.2f}} s")
'''
    scriptPath = jobDir / "reproDoranetJob.py"
    scriptPath.write_text(textwrap.dedent(scriptText), encoding="utf-8")
    scriptPath.chmod(0o755)


def runOneJob(
    jobIndex: int,
    starterSmiles: str,
    targetSmiles: str,
    generationRun: int,
    config: dict[str, Any],
    outputDir: Path,
) -> dict[str, Any]:
    startTime  = time.time()
    genSuffix  = str(config.get("generationSuffixTemplate", "_wGen{generation}")).format(generation=generationRun)
    jobName    = f"{config['baseJobPrefix']}{jobIndex}{genSuffix}"
    jobDir     = outputDir / jobName
    jobDir.mkdir(parents=True, exist_ok=True)

    atomMultiplier = float(config.get("atomMultiplier", 1.5))
    try:
        maxAtoms = computeMaxAtoms(starterSmiles, atomMultiplier)
    except ValueError as exc:
        return {
            "jobName": jobName, "jobDir": str(jobDir),
            "jobIndex": jobIndex, "starterSmiles": starterSmiles,
            "targetSmiles": targetSmiles, "generationRun": generationRun,
            "seconds": round(time.time() - startTime, 2),
            "status": "fail", "error": str(exc),
        }

    saveReproScript(jobDir, config, jobName, starterSmiles, targetSmiles, maxAtoms, generationRun)

    oldCwd = Path.cwd()
    os.chdir(jobDir)

    try:
        doranetPath = Path(config["doranetPath"]).expanduser().resolve()
        sys.path.insert(0, str(doranetPath))
        import doranet.modules.enzymatic as enzymatic

        starters   = {starterSmiles}
        helpers    = set(config["helpers"])
        userTarget = {targetSmiles}

        network = enzymatic.generate_network(
            job_name=jobName,
            starters=starters,
            gen=generationRun,
            max_atoms=maxAtoms,
            direction=config.get("direction", "forward"),
            targets=userTarget,
            ruleset=config["ruleset"],
        )

        smilesList = list(starters) + [m.uid for m in network.mols if m.uid not in starters]
        with Path(f"{jobName}_molecules.csv").open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["SMILES", "isStarter", "molFormula", "molWeight", "numHeavyAtoms"])
            for s in smilesList:
                writer.writerow([s, s in starters, *getSmilesProps(s)])

        targetCanonical = canonicalSmiles(targetSmiles)
        targetInNetwork = targetCanonical in {canonicalSmiles(s) for s in smilesList} if targetCanonical else False

        postInfo = maybeRunPostProcessing(network, jobName, userTarget, starters, helpers, config, generationRun)

        return {
            "jobName": jobName, "jobDir": str(jobDir),
            "jobIndex": jobIndex,
            "starterSmiles": starterSmiles,
            "targetSmiles": targetSmiles,
            "targetCanonical": targetCanonical,
            "generationRun": generationRun,
            "maxAtoms": str(maxAtoms),
            "targetInNetwork": targetInNetwork,
            "numMolecules": len(smilesList),
            "seconds": round(time.time() - startTime, 2),
            "status": "ok",
            **postInfo,
        }

    except Exception as exc:
        return {
            "jobName": jobName, "jobDir": str(jobDir),
            "jobIndex": jobIndex,
            "starterSmiles": starterSmiles,
            "targetSmiles": targetSmiles,
            "generationRun": generationRun,
            "maxAtoms": str(maxAtoms),
            "seconds": round(time.time() - startTime, 2),
            "status": "fail", "error": str(exc),
        }

    finally:
        os.chdir(oldCwd)


def validateConfig(config: dict[str, Any]) -> None:
    requiredKeys = ["doranetPath", "baseJobPrefix", "numParallelJobs", "helpers", "ruleset", "summaryCsv"]
    missing = [k for k in requiredKeys if k not in config]
    if missing:
        raise ValueError(f"Missing config keys: {missing}")
    if not resolveGenerationRuns(config):
        raise ValueError("No valid generations. Set generationsToRun, e.g. [1, 2, 3].")
    if not config.get("targetsCsvPath"):
        raise ValueError("targetsCsvPath must be set — starters and targets are read from the CSV.")


def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python doranet_pathways_multiGen.py config_doranet_multiGen.yaml")
        sys.exit(1)

    setThreadEnv()
    mainStart  = time.time()
    configPath = Path(sys.argv[1]).expanduser().resolve()
    config     = loadYaml(configPath)
    validateConfig(config)

    outputDir = Path(config.get("outputDir", ".")).expanduser().resolve()
    outputDir.mkdir(parents=True, exist_ok=True)

    csvPath      = Path(config["targetsCsvPath"]).expanduser().resolve()
    starterCol   = config.get("startersColumn", "starter_Canonical_SMILES")
    targetCol    = config.get("targetsColumn", "Canonical_SMILES")
    targetsToTest = config.get("targetsToTest", None)
    if targetsToTest is not None:
        targetsToTest = [int(x) for x in targetsToTest]

    pairs      = loadStarterTargetPairs(csvPath, starterCol, targetCol, targetsToTest)
    generations = resolveGenerationRuns(config)

    workItems = [
        (jobIndex, starterSmiles, targetSmiles, generationRun)
        for jobIndex, (starterSmiles, targetSmiles, _) in enumerate(pairs, start=1)
        for generationRun in generations
    ]

    atomMultiplier = float(config.get("atomMultiplier", 1.5))
    print(f"Config          : {configPath}")
    print(f"Output dir      : {outputDir}")
    print(f"CSV             : {csvPath}")
    print(f"Starter column  : {starterCol}")
    print(f"Target column   : {targetCol}")
    print(f"Unique pairs    : {len(pairs)}")
    print(f"Generations     : {generations}")
    print(f"Atom multiplier : {atomMultiplier}x (ceil)")
    print(f"Total jobs      : {len(workItems)} with nJobs={config['numParallelJobs']}")

    # Print maxAtoms preview per unique starter
    uniqueStarters = {s for s, _, _ in pairs}
    print("\nmaxAtoms preview per starter:")
    for smi in sorted(uniqueStarters):
        try:
            ma = computeMaxAtoms(smi, atomMultiplier)
            print(f"  {smi[:60]} → {ma}")
        except ValueError as exc:
            print(f"  {smi[:60]} → ERROR: {exc}")

    results = Parallel(n_jobs=int(config["numParallelJobs"]), backend="loky", verbose=10)(
        delayed(runOneJob)(jobIndex, starterSmiles, targetSmiles, generationRun, config, outputDir)
        for jobIndex, starterSmiles, targetSmiles, generationRun in workItems
    )

    resultsDf = pd.DataFrame(results)
    summaryPath = outputDir / config["summaryCsv"]
    resultsDf.to_csv(summaryPath, index=False)

    print(f"\nSaved summary : {summaryPath}")
    if "status" in resultsDf.columns:
        print(resultsDf["status"].value_counts(dropna=False).to_string())
    if {"generationRun", "status"}.issubset(resultsDf.columns):
        print("\nStatus by generation:")
        print(resultsDf.groupby(["generationRun", "status"]).size().unstack(fill_value=0).to_string())
    print(f"Total time    : {time.time() - mainStart:.2f} s")


if __name__ == "__main__":
    main()