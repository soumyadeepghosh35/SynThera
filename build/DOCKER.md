# Running SynThera with Docker

Everything needed to build and run SynThera in a container lives in this `build/`
directory. The image ships both conda environments and their Jupyter kernels
prebuilt, so a user only needs Docker. The rest of the repository is untouched;
it is bind-mounted into the container at runtime.

## How it works

The build context is `build/` only, so the image stays small and the large data
directories at the repository root never enter the context. The image contains:

1. JupyterLab in the base environment.
2. The `SynThera` (Steps 1 to 4) and `SynThera_ExptDesigns` (Step 5)
   environments, built from the pinned specs by the project's own `setup.py`.
3. Both kernels, registered so a single JupyterLab serves them.
4. RetroTide (JBEI Biosynthetic Cluster Simulator), installed into
   `SynThera_ExptDesigns`.

## Bundled tools versus ART

Every tool the pipeline needs is bundled in the image except ART. RetroTide is
open source, so it is installed into `SynThera_ExptDesigns` at build time with
`mapchiral` (its one otherwise-missing import) and no user action is required;
`import retrotide` and `import bcs` work on that kernel out of the box. Only ART
is licensed separately and left out (see below). To pin RetroTide for a
reproducible rebuild, pass `--build-arg RETROTIDE_REF="git+https://github.com/JBEI/RetroTide.git@<commit>"`.

At runtime, `docker-compose.yml` mounts the repository root at `/opt/synthera`,
so every notebook, `executables/` script, `external/ergochemics`, and any data
directory is live inside the container, and edits and outputs save back to the
host.

## Quick start

From this `build/` directory:

```bash
docker compose up --build
```

Watch the log for the tokenized URL (`http://127.0.0.1:8888/lab?token=...`) and
open it. In the launcher you will see both kernels, **SynThera** and
**SynThera_ExptDesigns**. Pick the one that matches the step you are running.
Everything except the ART-dependent step works at this point.

Stop with `docker compose down`.

## Three spec adjustments made for the build

`prepareEnvSpecs.py` rewrites copies of the specs inside the image, changing
exactly three pip entries and leaving every conda build pin untouched. Your
committed `environment.yml` and `exptDesigns_environment.yml` are not modified.

| Entry | Action | Reason |
|---|---|---|
| `art==2024.dev223+gbc27a92f1` | removed | ART is licensed separately (see below) and this pin is not on PyPI |
| `pathermo==0.0.1` | removed from both specs | not published to PyPI; the build cannot resolve it |
| `python-graphviz==0.21` | renamed to `graphviz==0.21` | the PyPI binding is named `graphviz`; `python-graphviz` is the conda name |

If a step needs `pathermo`, supply your own source the same way you would for ART.

## Installing ART

ART is free for non-commercial academic, non-profit, and U.S. government use but
licensed separately, so it is never baked into the image. `Step2_BuildModel`
needs it; every other step does not.

The full guide is in **ART.md**: how to obtain a license and the source, ART's
own maintainer-supported Docker image, and how to install ART into the SynThera
environment so `Step2_BuildModel` can import it. The short version of the
integrated route, once you are licensed:

```bash
# at build time, from build/
docker build \
  --build-arg ART_SOURCE="git+https://github.com/JBEI/AutomatedRecommendationTool.git" \
  -t synthera:with-art .

# or into a running container
./installArt.sh synthera "git+https://github.com/JBEI/AutomatedRecommendationTool.git"
docker commit synthera synthera:with-art
```

See ART.md for the licensing steps, the maintainer-supported route, and what to
do if a direct install does not resolve.

## Building by hand

```bash
docker build -t synthera:latest .

docker run --rm -p 8888:8888 -v "$(pwd)/..:/opt/synthera" synthera:latest
```

## Selecting the right kernel

| Notebook | Kernel |
|---|---|
| `Step1_DataPreparation.ipynb` | SynThera |
| `Step2_BuildModel.ipynb` | SynThera |
| `Step3_GenerateCompounds.ipynb` | SynThera |
| `Step4_PredictProperties.ipynb` | SynThera |
| `Step5_DesignExperiment.ipynb` | SynThera_ExptDesigns |

`Step0_GetBioactivityData.ipynb` runs on `SynThera`. Any ART-dependent step runs
only after ART is installed as described above.

## GPU and image size

The `SynThera` spec pins CUDA builds of `torch` and `tensorflow` with the full
NVIDIA runtime, so the image is large (on the order of ten-plus gigabytes) and a
GPU host with the NVIDIA container runtime is needed to use the GPU. On a CPU
host the packages install and fall back to CPU execution. To run on GPU, add
`--gpus all` to `docker run`, or uncomment the `deploy` block in
`docker-compose.yml`.

## Local (non-Docker) reproduction

To build the environments directly with conda instead of Docker, run the filter
first so the specs resolve, then the setup script:

```bash
cd build
python prepareEnvSpecs.py
python setup.py
```
