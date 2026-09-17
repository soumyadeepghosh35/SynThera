# Running SynThera with Docker

Everything needed to build and run SynThera in a container lives in this `build/`
directory. The image ships both conda environments and their Jupyter kernels
prebuilt, so a user only needs Docker. The rest of the repository is untouched;
it is bind-mounted into the container at runtime.

## How it works

The build context is `build/` only, so the image stays small and the large data
directories at the repository root never enter the context. The image contains:

1. JupyterLab in the base environment.
2. The `SynThera` (Steps 1, 2, and 4) and `SynThera_ExptDesigns` (Steps 3 and
   5) environments, built from the pinned specs by the project's own
   `setup.py`.
3. Both kernels, registered so a single JupyterLab serves them.
4. RetroTide (JBEI Biosynthetic Cluster Simulator), DORAnet (enzymatic
   reaction network generation), and DORA-XGB (enzymatic feasibility scoring),
   installed from their own GitHub repositories into whichever environment(s)
   use them.

## Bundled tools versus ART

Every tool the pipeline needs is bundled in the image except ART, and no user
action is required for any of them:

| Tool | Source | Installed into | Notes |
|---|---|---|---|
| RetroTide | [JBEI/RetroTide](https://github.com/JBEI/RetroTide) | `SynThera_ExptDesigns` | Installed with `mapchiral` (its one otherwise-missing import); `import retrotide` and `import bcs` work out of the box. |
| DORAnet | [wsprague-nu/doranet](https://github.com/wsprague-nu/doranet) | `SynThera` and `SynThera_ExptDesigns` | A real checkout is also kept on disk at `/opt/synthera-deps/doranet` (the `DORANET_DIR` build arg), because `runDORAnet.py` and the per-run pathway scripts require a directory at their configured `doranetPath`, not just an importable package. `executables/doranet_config.yaml` already points there. |
| DORA-XGB | [tyo-nu/DORA_XGB](https://github.com/tyo-nu/DORA_XGB) | `SynThera_ExptDesigns` only | Pins `xgboost==1.6.2`, which only this environment provides. |

Only ART is licensed separately and left out (see below). To pin any of these
to a specific commit for a reproducible rebuild, pass the matching build arg:

```bash
docker compose build \
  --build-arg RETROTIDE_REF="git+https://github.com/JBEI/RetroTide.git@<commit>" \
  --build-arg DORANET_REF="<commit-or-tag>" \
  --build-arg DORAXGB_REF="git+https://github.com/tyo-nu/DORA_XGB.git@<commit>"
```

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

## Spec adjustments made for the build

`prepareEnvSpecs.py` rewrites copies of the specs inside the image, leaving
every conda build pin untouched. Your committed `environment.yml` and
`exptDesigns_environment.yml` are not modified.

| Entry | Action | Reason |
|---|---|---|
| `art==2024.dev223+gbc27a92f1` | removed | ART is licensed separately (see below) and this pin is not on PyPI |
| `pathermo==0.0.1` | removed from both specs | not published to PyPI; the build cannot resolve it |
| `doranet==0.5.7a1` / `doranet==0.5.6a1` | removed from both specs | reinstalled afterward from its own GitHub source into both environments (see the table above), so the PyPI pin is never installed just to be overwritten |
| `dora-xgb==1.10` | removed from `exptDesigns_environment.yml` | reinstalled afterward from its own GitHub source, for the same reason |
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
| `Step3_GenerateCompounds.ipynb` | SynThera_ExptDesigns |
| `Step4_PredictProperties.ipynb` | SynThera |
| `Step5_DesignExperiment.ipynb` | SynThera_ExptDesigns |

`Step0_GetBioactivityData.ipynb` runs on `SynThera`. Any ART-dependent step runs
only after ART is installed as described above.

The `executables/` scripts need the same split: run `runADMET.py` and
`runART.py` with `conda run -n SynThera`, and `runDORAnet.py` and
`runDORAXGB.py` with `conda run -n SynThera_ExptDesigns` (the environment
Step 3 and Step 5 actually run on; it is also the only environment with the
`xgboost==1.6.2` that `runDORAXGB.py` needs). See the top-level README for the
full walkthrough.

## Apple Silicon (and any arm64 host)

Both environment specs are fully build-pinned **linux-64 (x86_64)** exports -
every package locked to an exact conda build hash. Docker builds against the
host's native architecture by default, so on an Apple Silicon Mac (arm64) the
base image resolves to `linux-aarch64` and the very first `mamba env create`
fails, since none of the pinned x86_64 build hashes exist in that channel. The
Dockerfile's `FROM --platform=linux/amd64 ...` line and `docker-compose.yml`'s
`platform: linux/amd64` both force the whole build and container to run under
x86_64 emulation instead, so this is handled automatically - no action needed
beyond having Docker Desktop's emulation support available (it ships with
QEMU-based emulation by default; enabling **Settings -> General -> Use Rosetta
for x86/amd64 emulation on Apple Silicon** speeds this up noticeably). Expect
the build, and the containerized packages, to run slower than they would
natively as a result.

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
