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
# ART (Automated Recommendation Tool) is intentionally NOT installed. ART is
# free for non-commercial academic, non-profit, and U.S. government use but
# licensed separately. Obtain a license and the source first (see ART.md), then
# add it with the ART_SOURCE build argument below or with installArt.sh against a
# running container. Only Step2_BuildModel needs ART.

# Pin this tag to a specific release on Docker Hub for reproducible rebuilds.
FROM condaforge/miniforge3:latest

# JupyterLab lives in the base environment and serves both kernels.
RUN mamba install -n base -y jupyterlab && mamba clean -a -y

# Build the environments in a location separate from the runtime mount point,
# so nothing here is shadowed when the project is mounted at /opt/synthera.
WORKDIR /opt/build

# Only the setup files are copied in. Relative paths resolve against this build/
# directory, which is the build context.
COPY environment.yml exptDesigns_environment.yml setup.py prepareEnvSpecs.py ./

# Rewrite the specs in place for a clean from-PyPI build: drop ART (licensed)
# and pathermo (absent from PyPI), and correct the graphviz pip name. The conda
# build pins are preserved byte for byte.
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
