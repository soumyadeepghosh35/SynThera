#!/usr/bin/env python3

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import sys
import time
import argparse

import numpy as np
import pandas as pd
import cloudpickle
import yaml

sys.path.append(".")
sys.path.append("/users/sghosh6/DTRA_project/MACAW/AutomatedRecommendationTool")

try:
    import warning_utils
    warning_utils.filter_end_user_warnings()
except ModuleNotFoundError:
    print("warning_utils not found; continuing without custom warning filtering.")

from art.core import RecommendationEngine
import art.utility as utils


def parseArgs():
    parser = argparse.ArgumentParser(
        description="Run ART workflow using a YAML configuration file."
    )
    parser.add_argument(
        "configFile",
        help="Path to YAML config file. Example: python runART.py artConfig.yaml"
    )
    return parser.parse_args()


def loadConfig(configFile):
    with open(configFile, "r") as f:
        config = yaml.safe_load(f)

    if config is None:
        raise ValueError("Config file is empty.")

    return config


def getRequiredConfig(config, keys):
    missingKeys = [key for key in keys if key not in config]
    if missingKeys:
        raise ValueError(f"Missing required config keys: {missingKeys}")


def buildPaths(config):
    runDir = os.getcwd()

    baseDataDir = config["baseDataDir"]
    modelBuildingDataSubdir = config.get("modelBuildingDataSubdir", ".")
    modelBuildingDataDir = os.path.join(baseDataDir, modelBuildingDataSubdir)

    outputDir = os.path.join(runDir, config["outputDir"])
    os.makedirs(outputDir, exist_ok=True)

    jobName = config["jobName"]

    defaultTrainingCsvName = f"{jobName}{config['trainingCsvSuffix']}"
    defaultArtReadyCsvName = f"{jobName}{config['artReadyCsvSuffix']}"

    inputTrainingCsv = config.get("inputTrainingCsv") or os.path.join(
        modelBuildingDataDir, defaultTrainingCsvName
    )

    artReadyCsv = config.get("artReadyCsv") or os.path.join(
        outputDir, defaultArtReadyCsvName
    )

    artModelFile = config.get("artModelFile") or os.path.join(
        outputDir, config["artModelFileName"]
    )

    pcaTransformerFile = os.path.join(outputDir, f"{jobName}_pca_transformer.cpkl")

    return {
        "runDir": runDir,
        "baseDataDir": baseDataDir,
        "modelBuildingDataDir": modelBuildingDataDir,
        "outputDir": outputDir,
        "inputTrainingCsv": inputTrainingCsv,
        "artReadyCsv": artReadyCsv,
        "artModelFile": artModelFile,
        "pcaTransformerFile": pcaTransformerFile,
    }


def loadTrainingData(config, paths):
    print("Loading training data...")
    DF = pd.read_csv(paths["inputTrainingCsv"])

    rowLimit = config.get("limitTrainingRows")
    if rowLimit is not None:
        DF = DF.iloc[:rowLimit].copy()
        print(f"Using first {len(DF)} rows from training data")

    return DF


def getInputAndResponseColumns(config, DF):
    inputColumnPrefix = config["inputFeatureColumnPrefix"]
    responseColumns = config["responseColumns"]

    inputColumns = [col for col in DF.columns if col.startswith(inputColumnPrefix)]

    if len(inputColumns) == 0:
        raise ValueError(
            f"No input feature columns found with prefix '{inputColumnPrefix}'."
        )

    missingResponseColumns = [col for col in responseColumns if col not in DF.columns]
    if missingResponseColumns:
        raise ValueError(
            f"Missing response columns in training data: {missingResponseColumns}"
        )

    print(f"Detected {len(inputColumns)} input feature columns")
    print(f"Response columns: {responseColumns}")

    return inputColumns, responseColumns


def applyPcaReduction(config, trainingDF, inputColumns, paths):
    """
    Optionally reduce input features via PCA before passing to ART.

    High-dimensional embeddings (e.g. 768-dim ChemBERTa) cause ART's
    cross-validation workers to exhaust per-process memory and get killed
    by the OS (SIGABRT-6).  Projecting to pcaComponents dimensions first
    eliminates that pressure while preserving nearly all variance in the
    embedding space that is predictive of pPotency.

    The fitted PCA transformer is saved alongside the ART model so it can
    be reapplied to new compound batches at inference time.
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    pcaComponents = config.get("pcaComponents")
    if pcaComponents is None:
        print("pcaComponents not set — skipping PCA reduction.")
        return trainingDF, inputColumns

    rawFeatures = trainingDF[inputColumns].to_numpy()
    nSamples, nFeatures = rawFeatures.shape
    nComponents = min(pcaComponents, nSamples, nFeatures)

    print(f"Applying PCA: {nFeatures} features → {nComponents} components "
          f"(pcaComponents={pcaComponents}, n_samples={nSamples})")

    scaler = StandardScaler()
    scaledFeatures = scaler.fit_transform(rawFeatures)

    pca = PCA(n_components=nComponents, random_state=config.get("randomSeed", 42))
    reducedFeatures = pca.fit_transform(scaledFeatures)

    explainedVar = np.sum(pca.explained_variance_ratio_) * 100
    print(f"PCA retained {explainedVar:.1f}% of total variance "
          f"across {nComponents} components.")

    pcaColumnNames = [f"ARTfeat_PC{i+1:04d}" for i in range(nComponents)]
    pcaDF = pd.DataFrame(reducedFeatures, columns=pcaColumnNames, index=trainingDF.index)

    nonFeatureCols = [c for c in trainingDF.columns if c not in inputColumns]
    reducedTrainingDF = pd.concat([trainingDF[nonFeatureCols], pcaDF], axis=1)

    pcaBundle = {"scaler": scaler, "pca": pca, "pcaColumnNames": pcaColumnNames}
    with open(paths["pcaTransformerFile"], "wb") as f:
        cloudpickle.dump(pcaBundle, f)
    print(f"Saved PCA transformer to: {paths['pcaTransformerFile']}")

    return reducedTrainingDF, pcaColumnNames


def saveArtReadyCsv(trainingDF, inputColumns, responseColumns, paths):
    print(f"Saving ART-ready CSV to: {paths['artReadyCsv']}")

    features = trainingDF[inputColumns].to_numpy()
    response = trainingDF[responseColumns].to_numpy()

    utils.save_edd_csv(
        features,
        response,
        inputColumns,
        paths["artReadyCsv"],
        responseColumns
    )

    print("Reloading ART-ready CSV...")
    artReadyDF = pd.read_csv(paths["artReadyCsv"])
    return artReadyDF


def loadArtModel(paths):
    print(f"Loading ART model from: {paths['artModelFile']}")
    with open(paths["artModelFile"], "rb") as f:
        artModel = cloudpickle.load(f)
    return artModel


def trainArtModel(config, artReadyDF, inputColumns, responseColumns, paths):
    # Cap CV worker count independently of artMaxMcmcCores.
    # joblib's loky backend sizes its pool from LOKY_MAX_CPU_COUNT when set.
    # Without this, CV forks as many workers as there are CPU cores, which
    # multiplies memory pressure by n_cores rather than by n_cv_folds.
    cvMaxWorkers = config.get("cvMaxWorkers", 4)
    os.environ["LOKY_MAX_CPU_COUNT"] = str(cvMaxWorkers)
    print(f"CV parallelism capped at {cvMaxWorkers} workers (cvMaxWorkers).")

    artParams = {
        "input_vars": inputColumns,
        "response_vars": responseColumns,
        "objective": config["artObjective"],
        "threshold": config["artThreshold"],
        "alpha": config["artAlpha"],
        "num_recommendations": config["artNumRecommendations"],
        "max_mcmc_cores": config["artMaxMcmcCores"],
        "seed": config["randomSeed"],
        "output_dir": paths["outputDir"],
        "recommend": config["artGenerateRecommendations"],
        "cross_val": config["artRunCrossValidation"],
        "num_tpot_models": config["artNumTpotModels"],
    }

    print("Training ART model...")
    startTime = time.time()
    artModel = RecommendationEngine(df=artReadyDF, **artParams)
    elapsedTime = time.time() - startTime
    print(f"ART training completed in {elapsedTime:.2f} seconds")

    if config["saveTrainedArtModel"]:
        with open(paths["artModelFile"], "wb") as f:
            cloudpickle.dump(artModel, f)
        print(f"Saved ART model to: {paths['artModelFile']}")

    return artModel


def getArtModel(config, artReadyDF, inputColumns, responseColumns, paths):
    if config["useExistingArtModel"]:
        return loadArtModel(paths)

    return trainArtModel(config, artReadyDF, inputColumns, responseColumns, paths)


def main():
    startTimeTotal = time.time()

    args = parseArgs()
    config = loadConfig(args.configFile)

    requiredKeys = [
        "baseDataDir",
        "modelBuildingDataSubdir",
        "outputDir",
        "jobName",
        "trainingCsvSuffix",
        "artReadyCsvSuffix",
        "artModelFileName",
        "inputFeatureColumnPrefix",
        "responseColumns",
        "useExistingArtModel",
        "saveTrainedArtModel",
        "artObjective",
        "artThreshold",
        "artAlpha",
        "artNumRecommendations",
        "artMaxMcmcCores",
        "artGenerateRecommendations",
        "artRunCrossValidation",
        "artNumTpotModels",
        "randomSeed",
    ]
    getRequiredConfig(config, requiredKeys)

    paths = buildPaths(config)

    print("Resolved paths:")
    for key, value in paths.items():
        print(f"  {key}: {value}")

    trainingDF = loadTrainingData(config, paths)
    inputColumns, responseColumns = getInputAndResponseColumns(config, trainingDF)

    # --- PCA reduction (new step) ---
    trainingDF, inputColumns = applyPcaReduction(config, trainingDF, inputColumns, paths)

    artReadyDF = saveArtReadyCsv(trainingDF, inputColumns, responseColumns, paths)
    artModel = getArtModel(config, artReadyDF, inputColumns, responseColumns, paths)

    elapsedTimeTotal = time.time() - startTimeTotal
    print(f"Total time: {elapsedTimeTotal:.2f} seconds")
    print(f"Job complete ...")

    return artModel


if __name__ == "__main__":
    main()
