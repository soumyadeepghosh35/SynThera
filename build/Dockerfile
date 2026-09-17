# SynThera reproducibility image
# Builds both conda environments and registers their Jupyter kernels so an end
# user only needs Docker, not a local conda toolchain. The build-pinned linux-64
# specs resolve as intended inside this linux-64 image.
#
# The build context is this build/ directory only, so the image stays small and
# the project's large data directories never enter the context. The project
# itself is bind-mounted at /opt/synthera at runtime (see docker-compose.yml),
# which keeps notebooks, executables, external/ergochemics, and any mounted data
# live inside the container with edits saved back to the host.
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
FROM condaforge/miniforge3:latest

# git is needed for the git+ pip installs below (RetroTide, and ART when
# licensed). Installing it at system level keeps it on PATH under any conda env.
RUN apt-get update && apt-get install -y --no-install-recommends git \
 && rm -rf /var/lib/apt/lists/*

# JupyterLab lives in the base environment and serves both kernels.
RUN mamba install -n base -y jupyterlab && mamba clean -a -y

# Build the environments in a location separate from the runtime mount point,
# so nothing here is shadowed when the project is mounted at /opt/synthera.
WORKDIR /opt/build

# Only the setup files are copied in. Relative paths resolve against this build/
# directory, which is the build context.
COPY environment.yml exptDesigns_environment.yml setup.py prepareEnvSpecs.py ./

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

# The project is bind-mounted here at runtime. Left empty in the image.
WORKDIR /opt/synthera

EXPOSE 8888

CMD ["jupyter", "lab", \
     "--ip=0.0.0.0", "--port=8888", \
     "--no-browser", "--allow-root", \
     "--ServerApp.root_dir=/opt/synthera"]
