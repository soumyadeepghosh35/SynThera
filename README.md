# SynThera

SynThera is an active-learning pipeline for AI-driven therapeutic design against
high-consequence pathogens. It couples machine-learned potency prediction,
retro-biosynthetic route generation, ADMET-based prioritization, thermodynamic
feasibility filtering, and host-optimized DNA design into a single
Design-Build-Test-Learn (DBTL) loop. The repository ships the notebooks,
reusable scripts, pinned environment specifications, and a Docker build so the
pipeline can be reproduced end to end, including by someone who has never used
conda before.

**New to this repository? Skip straight to [Quick start: Docker](#quick-start-docker-recommended).**
It is the supported way to run SynThera and the only path that does not require
installing or troubleshooting conda yourself.

## Pipeline at a glance

The workflow follows the DBTL framework and runs in five stages, executed once
per pathogen (case study):

1. **Step1_DataPreparation** - assemble and clean the bioactivity training data.
2. **Step2_BuildModel** - generate molecular embeddings and train/test potency
   predictors, including ART (Automated Recommendation Tool) models across
   several embedding methods.
3. **Step3_GenerateCompounds** - expand starting materials into novel candidates
   with DORAnet enzymatic retro-biosynthesis and RetroTide, then filter for
   drug-likeness.
4. **Step4_PredictProperties** - score the novel candidates for potency and
   ADMET properties.
5. **Step5_ExperimentalDesigns** - turn the prioritized candidates into
   experiment-ready, host-optimized DNA designs.

Steps 1-4 form the design-and-prioritization loop; Step 5 is the hand-off to
the wet lab.

## Repository layout

This is what actually exists in the repository today - not an idealized
template. Read this section before looking for a file by guessing its name.

| Path | What it is |
|---|---|
| `build/` | Everything needed to build and run the Docker image. See [Quick start](#quick-start-docker-recommended) and `build/DOCKER.md`. |
| `BacillusAnthracis/`, `Ebola/` | One completed case study per pathogen. Each contains its own `Step1_DataPreparation.ipynb` through `Step5_ExperimentalDesigns.ipynb`, and a `Results_<Pathogen>/` directory holding every output the run produced (embeddings, ART/ADMET predictions, pathway figures, DNA designs, and per-step subfolders with the config and log file used for that run). **Start here if you want to see or reproduce a finished result.** |
| `BioPKS_Compounds/` | Supplementary retro-biosynthesis exploration (RetroTide + DORAnet) for candidate antibiotic and antiviral scaffolds. This is prototyping work that feeds ideas into the pathogen case studies above; it is not itself a Step1-5 DBTL run and does not follow the same structure. |
| `Step0_GetBioactivityData.ipynb`, `Step1_DataPreparation.ipynb` | Shared, pathogen-agnostic data assembly, run once against the databases below. |
| `Step2_BuildModel.ipynb` ... `Step5_DesignExperiment.ipynb` (repository root) | **Empty template notebooks**, not finished runs. Copy one into a new pathogen folder as the starting point for a new case study; do not expect them to contain results. |
| `Input_Data/` | Shared reference data every case study reads: ChEMBL/PubChem exports, DrugBank references, DORAnet rulesets and enzyme-reaction caches, and the HuggingFace embedding models (ChemBERTa, MIST). Keep this alongside the notebooks. |
| `Database_Data_Raw/`, `Database_Data_Cleaned/` | The raw ChEMBL/PubChem pulls and their cleaned, per-pathogen bioactivity tables that Step0/Step1 produce and consume. |
| `ART_Results/` | Saved ART models per embedding method (ChemBERTa, MACAW, MIST, MolFormer), shared reference points for Step2. |
| `executables/` | Reusable, config-driven CLI scripts that the long-running steps call: `runDORAnet.py`, `runDORAXGB.py`, `runADMET.py`, `runART.py`, plus a `*_config.yaml` template for each. See [Running the long compute steps](#running-the-long-compute-steps). |
| `external/ergochemics/` | A bundled third-party cheminformatics library (reaction fingerprinting) that every Step5 notebook imports. See the caveat under [Known limitations](#known-limitations-read-before-you-build). |

## The two environments

The pipeline deliberately uses two separate conda environments so the heavier
design-and-modeling dependencies stay isolated from the DNA-design toolchain.
Both are built inside the Docker image automatically; you never need to create
them by hand unless you are doing a manual, non-Docker install.

| Environment | Spec file (in `build/`) | Covers | Key tooling |
|---|---|---|---|
| **SynThera** | `environment.yml` | Steps 1-4 | Embedding models, ART, DORAnet + RetroTide generation, potency and ADMET scoring |
| **SynThera_ExptDesigns** | `exptDesigns_environment.yml` | Step 5 | DNA Chisel, codon tables, Biopython, COBRApy, eQuilibrator, DORAnet, DORA-XGB, XGBoost, scikit-learn, RDKit |

Both spec files now declare their environment name internally (`name: SynThera`
/ `name: SynThera_ExptDesigns`), so this is what you get whether you build with
Docker, with `setup.py`, or by running `conda env create -f <file>` directly
with no extra flags.

## Quick start: Docker (recommended)

This is the path for a non-expert user. It needs Docker and nothing else - no
conda, no manually resolving Python packages.

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/)
   (macOS/Windows) or Docker Engine + Compose (Linux), and make sure it is
   running.
2. Open a terminal in this repository and go to the build directory:
   ```bash
   cd build
   docker compose up --build
   ```
   The first build downloads and compiles both environments and can take a
   while (the image is large, on the order of ten-plus gigabytes, because it
   pins CUDA builds of torch and tensorflow - see `build/DOCKER.md`).
3. Watch the terminal for a line like
   `http://127.0.0.1:8888/lab?token=...` and open that URL in your browser.
4. In the JupyterLab launcher you will see two kernels, **SynThera** and
   **SynThera_ExptDesigns**. Open a notebook and pick the kernel that matches
   the step (table below) from **Kernel -> Change Kernel** if it is not
   already selected.
5. Navigate to `BacillusAnthracis/` or `Ebola/` to open a finished case study,
   or copy the root-level template notebooks into a new folder to start one
   for a different pathogen.
6. Stop the container with `docker compose down` when you are done. Your
   notebooks, results, and edits are saved back to this folder on your
   computer the whole time - the container never holds the only copy.

Full details, GPU setup, and the ART installation steps below are in
[`build/DOCKER.md`](build/DOCKER.md).

### Selecting the right kernel

| Notebook | Kernel |
|---|---|
| `Step0_GetBioactivityData.ipynb` | SynThera |
| `Step1_DataPreparation.ipynb` | SynThera |
| `Step2_BuildModel.ipynb` | SynThera |
| `Step3_GenerateCompounds.ipynb` | SynThera |
| `Step4_PredictProperties.ipynb` | SynThera |
| `Step5_ExperimentalDesigns.ipynb` | SynThera_ExptDesigns |

### Installing ART (optional, only needed for Step2)

ART is free for non-commercial academic, non-profit, and U.S. government use,
but it is licensed separately from SynThera, so it is **not** baked into the
image by default. Every step except Step2's ART modeling works without it.
Once you have a license (see `build/ART.md`), install it with:

```bash
cd build
docker build --build-arg ART_SOURCE="git+https://github.com/JBEI/AutomatedRecommendationTool.git" -t synthera:with-art .
```

or into an already-running container with `build/installART.sh`. `build/ART.md`
has the full licensing steps and a fallback route if the direct install does
not resolve.

## Running the long compute steps

Some Step 3/4/5 work (large DORAnet pathway searches, ADMET scoring batches,
DORA-XGB feasibility ranking, ART cross-validation) is driven by standalone
scripts in `executables/`, not run directly inside a notebook cell, because it
is too long-running for an interactive kernel. The pattern is the same for
each:

1. The notebook that reaches that step writes out a small YAML config next to
   its results (for example `Results_<Pathogen>/Step4_ADMET/ADMET_config_....yaml`).
2. You run the matching script against that config, either from a terminal
   tab in JupyterLab (**File -> New -> Terminal**) or with
   `docker exec -it synthera bash`:
   ```bash
   conda run -n SynThera python executables/runDORAnet.py   <path-to-config>.yaml
   conda run -n SynThera python executables/runADMET.py     <path-to-config>.yaml
   conda run -n SynThera python executables/runDORAXGB.py   <path-to-config>.yaml
   conda run -n SynThera python executables/runART.py       <path-to-config>.yaml
   ```
3. The script prints progress and writes its results back into the same
   results folder; existing case studies keep the log as a `*.out` file next
   to the config for reference.

The `*_config.yaml` files directly under `executables/` are templates, not
live configs - copy one alongside a new run rather than editing it in place,
and see the caveat about hardcoded paths below before you rely on them.

## Known limitations (read before you build)

This section is intentionally blunt: these are gaps found by auditing the
repository, not yet fixed, and worth resolving (or at least being aware of)
before you rely on the container for a fresh reproduction.

- **Several `executables/*.yaml` templates hardcode the original author's
  filesystem paths.** For example, `executables/doranet_config.yaml` points
  `doranetPath` at `/users/sghosh6/DTRA_project/MACAW/doranet` and
  `executables/admet_config.yaml` points `inputFile`/`outputDir` at an NFS path
  under `/mnt/data.ese/...`. Neither path exists inside the container or on
  another machine. Edit these fields to match your own run before using a
  template config, and expect to do this per pathogen/target.
- **`external/ergochemics` is bundled but not installed as a package.** Step5
  notebooks add `external/ergochemics/src` to `sys.path` at runtime instead of
  a proper `pip install`. This works, but `ergochemics` declares
  `python>=3.12` and `rdkit>=2024.9.5` in its own `pyproject.toml`, while the
  `SynThera_ExptDesigns` environment it actually runs in pins Python 3.10.19
  and rdkit 2023.09.5. The notebooks already carry a compatibility monkeypatch
  around `Chem.MolToSmiles` for one API gap this causes; there is no guarantee
  it is the only one.
- **ART is excluded by design** (see above) - this is intentional, not a bug,
  but a non-expert user who opens `Step2_BuildModel.ipynb` without installing
  it first will hit an `ImportError` on `import art`.
- **The pinned pip stack in `environment.yml` is large and CUDA-specific**
  (torch, tensorflow, and a full NVIDIA CUDA 13 toolchain), and has only been
  validated on Linux x86-64. A clean `docker build` is the best way to confirm
  it still resolves; if a package fails to resolve, that is the first place to
  look.

## Manual (non-Docker) setup

If you would rather build the environments directly with conda:

```bash
cd build
python prepareEnvSpecs.py   # rewrites local copies of the specs for a from-PyPI build
python setup.py              # creates both environments and registers their kernels
```

`prepareEnvSpecs.py` removes the `art` and `pathermo` pip entries (neither
resolves from a clean PyPI install - `pathermo` is not currently imported by
any notebook or script in this repository, so removing it is safe) and renames
`python-graphviz` to the PyPI name `graphviz`. It edits copies of the specs in
place; if you run it directly against `build/environment.yml`, keep a backup
first.

Useful `setup.py` options:

```bash
python setup.py --force            # remove and recreate environments that already exist
python setup.py --skip-kernels     # build the environments without registering kernels
python setup.py --only SynThera    # act on a single environment
```

### Prerequisites for manual setup

- **conda** (Miniconda or Miniforge). [Mamba](https://mamba.readthedocs.io/) is
  used automatically when present, since it solves environments much faster.
- **Python 3** to run `setup.py` itself (the environments it creates pin their
  own Python versions - 3.12 for SynThera, 3.10 for SynThera_ExptDesigns).
- A **Linux x86-64** host for exact reproduction - both spec files are fully
  build-pinned linux-64 exports. On macOS or Windows, use the Docker route
  instead, or see the reproducibility note below.

### Reproducibility and platform notes

Both specification files are fully build-pinned linux-64 exports: every
package is fixed to an exact version and build hash, which is what makes them
reproduce bit-for-bit on Linux but also means the exact builds will not solve
on macOS or Windows. Options on a non-Linux machine:

- **Recommended:** use the Docker image above, which is already a Linux image.
- Recreate the environments from history on your own platform
  (`conda env export --from-history`) and let conda re-solve compatible
  builds. This trades exact reproduction for portability and drops the `pip:`
  section, so you will need to reattach the pip packages by hand.
- Rebuild from a loosely pinned spec (versions only, no build hashes) if you
  are prepared to validate that results match.

## Input data

`Input_Data/` holds everything the notebooks read to reproduce results:
ChEMBL/PubChem exports, DrugBank reference sets, DORAnet rulesets and
enzyme-reaction caches, and the HuggingFace embedding models. Keep it
alongside the notebooks so relative paths resolve. If it is distributed
separately because of size, download it and place it at the repository root
before running Step 1.

## Troubleshooting

- **`docker: command not found`** - install Docker Desktop (macOS/Windows) or
  Docker Engine (Linux) and make sure it is running before `docker compose up`.
- **`conda: command not found`** (manual setup only) - install Miniconda or
  Miniforge and reopen your shell so conda is on `PATH`.
- **The solve hangs or fails on macOS or Windows** (manual setup only) - this
  is the linux-64 pinning described above; use the Docker route or recreate
  from history.
- **`import art` fails** - expected until you install ART; see
  [Installing ART](#installing-art-optional-only-needed-for-step2).
- **A script under `executables/` cannot find its input file** - check the
  paths in the config YAML you passed it; several templates ship with the
  original author's local paths (see [Known limitations](#known-limitations-read-before-you-build)).
- **The kernel doesn't appear in Jupyter** (manual setup only) - confirm the
  environment was created, then rerun `python setup.py --only <env>` or the
  manual `ipykernel install` command from `build/DOCKER.md`.
- **A pip package fails to build** - check that system build tools are
  available; on a fresh Linux host you may need `build-essential`. The Docker
  route avoids this entirely.
