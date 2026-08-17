# SynThera

SynThera is an end-to-end, active-learning pipeline for AI-driven therapeutic design against high-consequence pathogens. It couples machine-learned potency prediction, retro-biosynthetic route generation, ADMET-based prioritization, thermodynamic feasibility filtering, and host-optimized DNA design into a single Design–Build–Test–Learn (DBTL) loop. The repository ships the notebooks, scripts, and pinned environment specifications needed to reproduce the computational results end to end.

This README explains the pipeline, the two computing environments it needs, and how to set everything up with a single command.

## Pipeline at a glance

The workflow follows the DBTL framework and runs in five stages. Steps 1 through 4 form the **design and prioritization** loop; Step 5 turns the surviving candidates into **experiment-ready DNA constructs**.

1. **Data preparation** — assemble and clean the training data for the potency models.
2. **Model building** — generate molecular embeddings and train and test the potency predictors, including the ART (Automated Recommendation Tool) models across several embedding methods.
3. **Compound generation** — expand starting materials into novel candidates using DORAnet enzymatic retro-biosynthesis and RetroTide, then filter for drug-likeness.
4. **Property prediction** — score the novel candidates for potency and ADMET properties.
5. **Experiment design** — deliver host-optimized, experiment-ready DNA designs for the prioritized candidates.

## Repository structure

The table below mirrors the code structure of the project (Supplementary Table S1).

| File | Type | Purpose |
|---|---|---|
| `environment.yml` | Environment | Reproduces the `SynThera` environment for Steps 1–4 |
| `exptDesigns_environment.yml` | Environment | Reproduces the `SynThera_ExptDesigns` environment for Step 5 |
| `Step1_DataPreparation.ipynb` | Notebook | Prepares the data used to train the models |
| `Step2A_BuildModel.ipynb` | Notebook | Generates embeddings, then trains and tests the models |
| `Step2B_BuildModel_wART` | Python | Reproduces the ART models across the different embeddings |
| `Step3A_GenerateCompounds` | Python | Generates candidates via DORAnet + RetroTide |
| `Step3B_GenerateCompounds.ipynb` | Notebook | Filters the novel, drug-like compounds from retro-biosynthesis |
| `Step4_PredictProperties.ipynb` | Notebook | Predicts potency and ADMET properties of the novel candidates |
| `Step5_DesignExperiment.ipynb` | Notebook | Delivers experiment-ready DNA designs |
| `Input_Data` | Directory | Input data required to reproduce the results |
| `setup.py` | Script | Builds both environments and registers their Jupyter kernels |

## The two environments

The pipeline deliberately uses two separate conda environments so that the heavier design-and-modeling dependencies stay isolated from the DNA-design toolchain, which pins its own versions of thermodynamics and sequence-engineering libraries.

| Environment | Covers | Key tooling |
|---|---|---|
| **SynThera** | Steps 1–4 | Embedding, ART modeling, DORAnet + RetroTide generation, potency and ADMET scoring |
| **SynThera_ExptDesigns** | Step 5 | DNA Chisel, codon tables, Biopython, COBRApy, eQuilibrator, DORAnet, DORA-XGB, XGBoost, scikit-learn, RDKit |

Each environment gets its own Jupyter kernel of the same name, so you can pick the right interpreter directly from the notebook.

## Prerequisites

- **conda** (Miniconda or Miniforge). [Mamba](https://mamba.readthedocs.io/) is recommended and, if present, is used automatically because it solves environments far faster.
- **Python 3** to run `setup.py` (any interpreter that can launch the script; the environments themselves pin Python 3.10).
- A **Linux x86-64** host for exact reproduction. See the reproducibility note below before running on macOS or Windows.

## Quick start

Place `setup.py`, `environment.yml`, and `exptDesigns_environment.yml` in the same directory, then run:

```bash
python setup.py
```

That single command creates both environments and registers both Jupyter kernels. When it finishes, launch Jupyter and select the kernel that matches the step you're running:

```bash
jupyter lab
```

Useful options:

```bash
python setup.py --force            # remove and recreate environments that already exist
python setup.py --skip-kernels     # build the environments without registering kernels
python setup.py --only SynThera    # act on a single environment
```

## What `setup.py` does

For each environment it:

1. Locates the solver, preferring `mamba` over `conda` when both are available.
2. Checks whether the environment already exists, and either leaves it in place or recreates it when `--force` is set.
3. Creates the environment from its specification file, overriding the internal `name:` field so the environments are named `SynThera` and `SynThera_ExptDesigns` regardless of what the exported file was called.
4. Ensures `ipykernel` is present, then registers a Jupyter kernel whose name and display label match the environment.

If a specification file is missing, the script skips that environment with a clear message rather than failing the whole run, so you can set up one environment at a time if needed.

## Selecting the right kernel

| Notebook / script | Kernel |
|---|---|
| `Step1_DataPreparation.ipynb` | SynThera |
| `Step2A_BuildModel.ipynb` | SynThera |
| `Step2B_BuildModel_wART` | SynThera |
| `Step3A_GenerateCompounds` | SynThera |
| `Step3B_GenerateCompounds.ipynb` | SynThera |
| `Step4_PredictProperties.ipynb` | SynThera |
| `Step5_DesignExperiment.ipynb` | SynThera_ExptDesigns |

In Jupyter, the kernel appears in the launcher and under **Kernel → Change Kernel** with the environment's name.

## Reproducibility and platform notes

Both specification files are **fully build-pinned linux-64 exports**: every package is fixed to an exact version and build hash. That's what makes them reproduce bit-for-bit on Linux, but it also means the exact builds won't solve on macOS or Windows. You have three options on a non-Linux machine:

- **Recommended:** run inside a Linux container or a Linux HPC node, where the pinned specs resolve as intended.
- Recreate the environments from history on your own platform (`conda env export --from-history`) and let conda re-solve compatible builds. This trades exact reproduction for portability, and it drops the `pip:` section, so you'll need to reattach the pip packages by hand.
- Rebuild from a loosely pinned spec (versions only, no build hashes) if you're prepared to validate that results match.

`setup.py` prints a warning when it detects a non-Linux host so you aren't caught out by a failed solve.

## Input data

The `Input_Data` directory holds everything the notebooks read to reproduce the results. Keep it alongside the notebooks so the relative paths resolve. If the directory is distributed separately because of size, download it and place it at the repository root before running Step 1.

## Manual setup without `setup.py`

If you'd rather set the environments up by hand:

```bash
# Steps 1-4
conda env create -n SynThera -f environment.yml
conda run -n SynThera python -m ipykernel install --user \
  --name SynThera --display-name SynThera

# Step 5
conda env create -n SynThera_ExptDesigns -f exptDesigns_environment.yml
conda run -n SynThera_ExptDesigns python -m ipykernel install --user \
  --name SynThera_ExptDesigns --display-name SynThera_ExptDesigns
```

## Troubleshooting

- **`conda: command not found`** — install Miniconda or Miniforge and reopen your shell so conda is on `PATH`.
- **The solve hangs or fails on macOS or Windows** — this is the linux-64 pinning described above; use a Linux environment or recreate from history.
- **The kernel doesn't appear in Jupyter** — confirm the environment was created, then rerun `python setup.py --only <env>` or the manual `ipykernel install` command. Kernels install per user under your Jupyter data directory.
- **A pip package fails to build during environment creation** — check that system build tools are available; on a fresh Linux host you may need `build-essential`. The container route avoids this entirely.
