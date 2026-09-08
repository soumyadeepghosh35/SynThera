#!/usr/bin/env python3
import os, time, argparse, yaml
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

os.environ.update({
    "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
})


def loadConfig(path):
    with open(path, "r") as f:
        cfg = yaml.safe_load(f)
    if cfg is None:
        raise ValueError("Config file is empty.")
    return cfg


def validateConfig(cfg):
    for k in ("inputFile", "SMILES_Column", "batchSize"):
        if k not in cfg:
            raise KeyError(f"Missing required config key: {k}")


def buildPaths(cfg):
    outputDir = Path(cfg.get("outputDir", ".")).resolve()
    inputPath = Path(cfg["inputFile"])
    if not inputPath.is_absolute():
        inputPath = (Path.cwd() / inputPath).resolve()
    if not inputPath.exists():
        raise FileNotFoundError(f"Input file not found: {inputPath}")

    stem = inputPath.stem
    outSuffix = cfg.get("outputSuffix", "_wADMET")
    divSuffix = cfg.get("diversityOutputSuffix", "_diverseInput")

    return {
        "inputCsvPath": str(inputPath),
        "outputCsvPath": str(outputDir / f"{stem}{outSuffix}.csv"),
        "checkpointCsvPath": str(outputDir / f"{stem}{outSuffix}_checkpoint.csv"),
        "diversityCsvPath": str(outputDir / f"{stem}{divSuffix}.csv"),
        "inputFile": inputPath.name,
    }


def loadAdmetModel():
    from admet_ai import ADMETModel
    return ADMETModel()


def readInputCsv(inputCsvPath, selectSMILES):
    df = pd.read_csv(inputCsvPath)
    if selectSMILES is not None:
        if not isinstance(selectSMILES, int) or selectSMILES <= 0:
            raise ValueError("selectSMILES must be a positive integer or null/None.")
        df = df.head(selectSMILES).copy()
    return df.reset_index(drop=True)


def canonicalizeSmiles(s):
    if pd.isna(s):
        return None
    mol = Chem.MolFromSmiles(s)
    return None if mol is None else Chem.MolToSmiles(mol, canonical=True)


def deduplicateForPrediction(df, smiles_col):
    return df[[smiles_col]].dropna().drop_duplicates().reset_index(drop=True)


def applyDiversityFilter(df, smiles_col, similarityThreshold, fingerprintRadius, fingerprintBits, useCanonicalFiltering, diversityChunkSize):
    """Greedy Tanimoto diversity filter; returns filtered full dataframe (all original columns)."""
    uniq = df[[smiles_col]].dropna().drop_duplicates().reset_index(drop=True)

    selectedSmiles, selectedCanonicalSmiles, selectedFingerprints = [], [], []
    seen = set()
    validMolCount = 0

    for start in range(0, len(uniq), diversityChunkSize):
        chunk = uniq.iloc[start:start + diversityChunkSize]
        before = len(selectedSmiles)

        for smi in chunk[smiles_col]:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            validMolCount += 1
            can = Chem.MolToSmiles(mol, canonical=True)
            if can in seen:
                continue

            fp = AllChem.GetMorganFingerprintAsBitVect(mol, fingerprintRadius, nBits=fingerprintBits)

            if not selectedFingerprints:
                selectedSmiles.append(smi)
                selectedCanonicalSmiles.append(can)
                selectedFingerprints.append(fp)
                seen.add(can)
            else:
                sims = DataStructs.BulkTanimotoSimilarity(fp, selectedFingerprints)
                if max(sims) < similarityThreshold:
                    selectedSmiles.append(smi)
                    selectedCanonicalSmiles.append(can)
                    selectedFingerprints.append(fp)
                    seen.add(can)

        after = len(selectedSmiles)
        print(f"diversity: processed {min(start + diversityChunkSize, len(uniq))}/{len(uniq)} | +{after - before} new | total {after}")

    diverseFiltered = pd.DataFrame({smiles_col: selectedSmiles, "CanonicalSMILES": selectedCanonicalSmiles})

    if useCanonicalFiltering:
        selSet = set(selectedCanonicalSmiles)
        mask = df[smiles_col].apply(canonicalizeSmiles).isin(selSet)
    else:
        selSet = set(selectedSmiles)
        mask = df[smiles_col].isin(selSet)

    filtered = df.loc[mask].reset_index(drop=True)
    stats = {
        "originalUniqueValidMolecules": validMolCount,
        "diversityFilteredMolecules": len(diverseFiltered),
        "finalFilteredRows": len(filtered),
    }
    return filtered, diverseFiltered, stats


def splitIntoBatches(items, batchSize):
    for i in range(0, len(items), batchSize):
        yield items[i:i + batchSize]


def safePredictBatch(model, smilesBatch, smiles_col):
    pred = model.predict(smiles=smilesBatch)
    if not isinstance(pred, pd.DataFrame):
        pred = pd.DataFrame(pred)
    pred = pred.reset_index(drop=True)
    if len(pred) != len(smilesBatch):
        raise ValueError("Prediction batch size mismatch.")
    return pd.concat([pd.DataFrame({smiles_col: smilesBatch}), pred], axis=1)


def loadCheckpoint(path, smiles_col):
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    # tolerate older checkpoints that used 'SMILES'
    if smiles_col not in df.columns and "SMILES" in df.columns:
        df = df.rename(columns={"SMILES": smiles_col})
    if smiles_col not in df.columns:
        raise ValueError(f"Checkpoint must contain '{smiles_col}' (or 'SMILES').")
    return df.drop_duplicates(subset=smiles_col, keep="last").reset_index(drop=True)


def writeBatchToCheckpoint(batchDF, path, writeHeader):
    Path(os.path.dirname(path)).mkdir(parents=True, exist_ok=True)
    batchDF.to_csv(path, mode="a", header=writeHeader, index=False)


def compactCheckpoint(path, smiles_col):
    df = pd.read_csv(path)
    if smiles_col not in df.columns and "SMILES" in df.columns:
        df = df.rename(columns={"SMILES": smiles_col})
    df = df.drop_duplicates(subset=smiles_col, keep="last").reset_index(drop=True)
    df.to_csv(path, index=False)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("configFile")
    args = ap.parse_args()

    t0 = time.time()
    cfg = loadConfig(args.configFile)
    validateConfig(cfg)
    paths = buildPaths(cfg)

    smiles_col = cfg["SMILES_Column"]
    if smiles_col is None or smiles_col == "":
        raise ValueError("SMILES_Column must be a non-empty string.")

    selectSMILES = cfg.get("selectSMILES", None)
    batchSize = int(cfg["batchSize"])

    saveEveryBatch = bool(cfg.get("saveEveryBatch", True))
    resumeFromCheckpoint = bool(cfg.get("resumeFromCheckpoint", True))
    deduplicateSmiles = bool(cfg.get("deduplicateSmiles", True))

    enableDiversityFilter = bool(cfg.get("enableDiversityFilter", False))
    similarityThreshold = float(cfg.get("similarityThreshold", 0.75))
    fingerprintRadius = int(cfg.get("fingerprintRadius", 2))
    fingerprintBits = int(cfg.get("fingerprintBits", 2048))
    diversityChunkSize = int(cfg.get("diversityChunkSize", 5000))
    saveDiversityFilteredInput = bool(cfg.get("saveDiversityFilteredInput", True))
    useCanonicalFiltering = bool(cfg.get("useCanonicalFiltering", True))

    print("\n=== ADMET profile prediction ===")
    print(f"inputFile: {paths['inputFile']}")
    print(f"SMILES column: {smiles_col}")
    print(f"outputCsv: {paths['outputCsvPath']}")
    print(f"checkpoint: {paths['checkpointCsvPath']}\n")

    df = readInputCsv(paths["inputCsvPath"], selectSMILES)
    if smiles_col not in df.columns:
        raise ValueError(f"Input CSV does not contain column '{smiles_col}'. Available columns: {list(df.columns)}")

    if enableDiversityFilter:
        print("Applying diversity filter...")
        df, diverseDF, stats = applyDiversityFilter(
            df, smiles_col,
            similarityThreshold, fingerprintRadius, fingerprintBits,
            useCanonicalFiltering, diversityChunkSize
        )
        print(f"diversity stats: {stats}")
        if saveDiversityFilteredInput:
            Path(os.path.dirname(paths["diversityCsvPath"])).mkdir(parents=True, exist_ok=True)
            df.to_csv(paths["diversityCsvPath"], index=False)

    predInput = deduplicateForPrediction(df, smiles_col) if deduplicateSmiles else df[[smiles_col]].dropna().reset_index(drop=True)
    allSmiles = predInput[smiles_col].tolist()
    print(f"SMILES for ADMET prediction: {len(allSmiles)}")

    predDF = None
    processed = set()

    if resumeFromCheckpoint:
        predDF = loadCheckpoint(paths["checkpointCsvPath"], smiles_col)
        if predDF is not None:
            processed = set(predDF[smiles_col].dropna().tolist())
            print(f"Loaded checkpoint: {len(predDF)} SMILES")

    remaining = [s for s in allSmiles if s not in processed]
    print(f"Remaining SMILES: {len(remaining)}")

    if remaining:
        model = loadAdmetModel()
        writeHeader = not os.path.exists(paths["checkpointCsvPath"]) or os.path.getsize(paths["checkpointCsvPath"]) == 0

        for i, batch in enumerate(splitIntoBatches(remaining, batchSize), start=1):
            t_batch = time.time()
            batchPred = safePredictBatch(model, batch, smiles_col)
            if saveEveryBatch:
                writeBatchToCheckpoint(batchPred, paths["checkpointCsvPath"], writeHeader)
                writeHeader = False
            else:
                raise RuntimeError("saveEveryBatch=False is not supported.")
            print(f"batch {i}: size={len(batch)} time={time.time() - t_batch:.2f}s")

        predDF = compactCheckpoint(paths["checkpointCsvPath"], smiles_col)
    else:
        if predDF is None:
            raise RuntimeError("No SMILES processed and no checkpoint found.")
        predDF = predDF.copy()

    # Merge predictions back; then append prediction columns to the end
    orig_cols = df.columns.tolist()
    pred_cols = [c for c in predDF.columns if c != smiles_col]

    merged = df.merge(predDF, on=smiles_col, how="left")
    final_order = orig_cols + [c for c in pred_cols if c in merged.columns]
    # keep any unexpected columns at the end too
    final_order += [c for c in merged.columns if c not in final_order]
    merged = merged[final_order]

    Path(os.path.dirname(paths["outputCsvPath"])).mkdir(parents=True, exist_ok=True)
    merged.to_csv(paths["outputCsvPath"], index=False)

    print(f"\nSaved: {paths['outputCsvPath']}")
    print(f"Total runtime: {time.time() - t0:.2f}s")


if __name__ == "__main__":
    main()