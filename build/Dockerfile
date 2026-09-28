# SynThera reproducibility image
# Builds both conda environments and registers their Jupyter kernels so an end
# user only needs Docker, not a local conda toolchain. The build-pinned linux-64
# specs resolve as intended inside this linux-64 image.
#
# The build context is the repository root (see docker-compose.yml's
# build.context), filtered through the root ../.dockerignore, which excludes
# everything except build/, executables/, tests/, and external/ - so the
# project's multi-GB per-pathogen data directories never enter the context or
# the image. executables/, tests/, and external/ergochemics/src are copied
# into a fixed, standalone location (/opt/synthera-app, see the bottom of this
# file) so that a user who only has the built/pulled image - no source
# checkout, no bind mount - can still run the test suite. At runtime the full
# project is additionally bind-mounted at /opt/synthera (see
# docker-compose.yml), which keeps notebooks, executables, external/ergochemics,
# and any mounted data live inside the container with edits saved back to the
# host; that live copy shadows the baked-in one for anyone using compose.
#
# RetroTide (JBEI Biosynthetic Cluster Simulator), DORAnet (enzymatic reaction
# network generation), and DORA-XGB (enzymatic feasibility scoring) are all
# open source and are baked into the image from their own GitHub repositories,
# so users install nothing for them. ART (Automated Recommendation Tool) is the
# one exception: it is licensed separately and is NOT installed. Obtain an ART
# license and source first (see ART.md), then add it with the ART_SOURCE build
# argument below or with installArt.sh against a running container. Only
# Step2_BuildModel needs ART.

# Pin this tag to a specific release on Docker Hub for reproducible rebuilds.
# --platform=linux/amd64 is required even on an Apple Silicon (arm64) host:
# environment.yml and exptDesigns_environment.yml are fully build-pinned
# linux-64 (x86_64) exports (exact conda build hashes), which do not exist in
# the linux-aarch64 channel. Without this, Docker pulls the base image for the
# host's native architecture and the very first `mamba env create` fails with
# "does not exist (perhaps a typo or a missing channel)" for nearly every
# package. This line forces the whole image to build and run under x86_64
# emulation instead, so the pinned specs resolve exactly as exported.
FROM --platform=linux/amd64 condaforge/miniforge3:latest

# Make every RUN below fail immediately if any command in a &&/; chain or a
# pipe fails, instead of only failing on the exit status of the very last
# command executed. Without this, a `for env in ...; do cmdA && cmdB; done`
# RUN step can report success even though cmdB genuinely failed for one of
# the environments - which is exactly what happened with the DORAnet import
# check below: it failed for SynThera, the loop carried on to
# SynThera_ExptDesigns anyway, that one succeeded, and the RUN exited 0.
SHELL ["/bin/bash", "-eo", "pipefail", "-c"]

# git is needed for the git+ pip installs below (RetroTide, DORAnet,
# DORA-XGB, and ART when licensed); installing it at system level keeps it
# on PATH under any conda env. libxrender1/libxext6/libsm6 are X11 shared
# libraries that rdkit's Chem.Draw submodule links against and dlopens at
# import time even when it never opens a display or draws anything - DORAnet
# imports Chem.Draw unconditionally as part of its own post_processing
# module, so without these three packages the import fails with
# "ImportError: libXrender.so.1: cannot open shared object file", which is
# the failure the SHELL fix above now actually surfaces instead of hiding.
RUN apt-get update && apt-get install -y --no-install-recommends \
      git \
      libxrender1 \
      libxext6 \
      libsm6 \
 && rm -rf /var/lib/apt/lists/*

# JupyterLab lives in the base environment and serves both kernels.
RUN mamba install -n base -y jupyterlab && mamba clean -a -y

# Build the environments in a location separate from the runtime mount point,
# so nothing here is shadowed when the project is mounted at /opt/synthera.
WORKDIR /opt/build

# Only the setup files are copied in here (executables/, tests/, and
# external/ are copied separately near the end of this file). Source paths
# are prefixed with build/ because the context is now the repository root.
COPY build/environment.yml build/exptDesigns_environment.yml build/setup.py build/prepareEnvSpecs.py ./

# Rewrite the specs in place for a clean from-PyPI build: drop ART (licensed),
# pathermo (absent from PyPI), and doranet/dora-xgb (installed below from their
# own GitHub sources instead, so the exact code the executables/ scripts expect
# is guaranteed), and correct the graphviz pip name. The conda build pins are
# preserved byte for byte.
RUN python prepareEnvSpecs.py

# Build both environments with the project's own logic. The kernel step is
# skipped here because it registers per user; kernels are registered into the
# base prefix below so they are visible to any user of the image.
RUN python setup.py --skip-kernels && mamba clean -a -y

# Ensure ipykernel is present in each environment, then register both kernels
# into the base prefix so the base JupyterLab discovers them.
RUN for env in SynThera SynThera_ExptDesigns; do \
      conda run -n "$env" python -c "import ipykernel" 2>/dev/null \
        || mamba install -n "$env" -y ipykernel; \
      conda run -n "$env" python -m ipykernel install \
        --prefix /opt/conda --name "$env" --display-name "$env"; \
    done && mamba clean -a -y

# pytest is test tooling, not a pipeline dependency, so it is deliberately not
# part of the pinned environment.yml / exptDesigns_environment.yml exports.
# It is installed into both environments here so tests/unit/ can run with
# `conda run -n <env> pytest tests/unit` without touching the pinned specs.
RUN for env in SynThera SynThera_ExptDesigns; do \
      conda run -n "$env" python -c "import pytest" 2>/dev/null \
        || mamba install -n "$env" -y pytest; \
    done && mamba clean -a -y

# RetroTide is installed into SynThera_ExptDesigns, which already provides its
# dependencies (cobra, rdkit, numpy, scipy, typing_extensions). mapchiral is a
# required top-level import that ships no dependencies of its own. Both are
# installed with --no-deps so the build-pinned environment is left untouched,
# and the import is verified so a broken install fails the build early. Pin
# RETROTIDE_REF to a commit (append @<sha>) for reproducible rebuilds.
ARG RETROTIDE_ENV="SynThera_ExptDesigns"
ARG RETROTIDE_REF="git+https://github.com/JBEI/RetroTide.git"
RUN conda run -n "$RETROTIDE_ENV" pip install --no-deps mapchiral \
 && conda run -n "$RETROTIDE_ENV" pip install --no-deps "$RETROTIDE_REF" \
 && conda run -n "$RETROTIDE_ENV" python -c "import retrotide, bcs; print('RetroTide import OK')"

# DORAnet (enzymatic reaction network generation) is what SynThera_ExptDesigns
# actually uses, for both Step 3 (compound generation - Step3_GenerateCompounds.ipynb
# runs on SynThera_ExptDesigns, not SynThera) and Step 5 (pathway design). The
# executables/ scripts that drive it (runDORAnet.py, and the per-run pathway
# scripts Step 5 generates) require a real on-disk checkout at the
# `doranetPath` set in their config, not just an importable package, so one
# checkout is kept in the image at DORANET_DIR. It is installed into SynThera
# too, out of caution in case anything in Steps 1/2/4 still depends on it -
# drop that half of the loop below if you confirm nothing does. Pass
# --build-arg DORANET_REF=<tag-or-commit> to pin a specific release instead of
# the default branch.
ARG DORANET_REPO="https://github.com/wsprague-nu/doranet.git"
ARG DORANET_REF=""
ARG DORANET_DIR="/opt/synthera-deps/doranet"
RUN git clone "$DORANET_REPO" "$DORANET_DIR" \
 && if [ -n "$DORANET_REF" ]; then git -C "$DORANET_DIR" checkout "$DORANET_REF"; fi \
 && for env in SynThera SynThera_ExptDesigns; do \
      conda run -n "$env" pip install --no-deps "$DORANET_DIR" \
       && conda run -n "$env" python -c "import doranet.modules.enzymatic, doranet.modules.post_processing; print('DORAnet import OK in $env')"; \
    done

# DORA-XGB (gradient-boosted enzymatic-reaction feasibility scoring) pins
# xgboost==1.6.2, which only SynThera_ExptDesigns provides, so it is installed
# there only. Pass --build-arg DORAXGB_REF=<git-spec-with-@commit> to pin a
# specific release instead of the default branch.
ARG DORAXGB_REF="git+https://github.com/tyo-nu/DORA_XGB.git"
RUN conda run -n SynThera_ExptDesigns pip install --no-deps "$DORAXGB_REF" \
 && conda run -n SynThera_ExptDesigns python -c "from DORA_XGB import DORA_XGB; print('DORA-XGB import OK')"

# Optional ART install once you are licensed. ART_SOURCE may be a git URL, a
# PyPI-style spec, or a wheel reachable during the build, for example:
#   git+https://github.com/JBEI/AutomatedRecommendationTool.git
# Empty by default, so the standard build ships without ART. See ART.md.
ARG ART_SOURCE=""
RUN if [ -n "$ART_SOURCE" ]; then \
      echo "Installing ART into SynThera from: $ART_SOURCE"; \
      conda run -n SynThera pip install "$ART_SOURCE"; \
    else \
      echo "ART not installed. Provide --build-arg ART_SOURCE=... once licensed."; \
    fi

# Bake a standalone copy of the small, non-data parts of the project into a
# fixed location so the test suite (and the scripts it exercises) works from
# just the image - no bind mount, no source checkout required. This is what
# `docker run --rm <image> bash /opt/synthera-app/tests/run_all_tests.sh`
# (see tests/README.md) runs against. Anyone using docker-compose still gets
# the live project bind-mounted at /opt/synthera, which shadows this copy.
COPY executables /opt/synthera-app/executables
COPY tests /opt/synthera-app/tests
COPY external /opt/synthera-app/external

# The project is bind-mounted here at runtime. Left empty in the image.
WORKDIR /opt/synthera

EXPOSE 8888

CMD ["jupyter", "lab", \
     "--ip=0.0.0.0", "--port=8888", \
     "--no-browser", "--allow-root", \
     "--ServerApp.root_dir=/opt/synthera"]
