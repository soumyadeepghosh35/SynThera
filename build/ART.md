# Getting and installing ART

SynThera's model-building step (`Step2_BuildModel.ipynb`) uses ART, the Automated
Recommendation Tool from Lawrence Berkeley National Laboratory. ART is **not**
bundled in the SynThera image, because it is licensed separately. This guide
explains how to obtain ART, the two supported ways to install it, and how it
plugs into the SynThera workflow.

Reference: Radivojevic T., Costello Z., Workman K., Garcia Martin H., "A machine
learning Automated Recommendation Tool for synthetic biology", Nature
Communications 11, 4879 (2020). Docs: https://lbl-biosci.gitlab.io/ese/art/

## 1. Get a license and the source

ART is free for non-commercial use by academic, non-profit, and U.S. government
institutions, and is distributed patent-pending. Before installing it, obtain a
license and access to the source code through the JBEI/ART project:

- License and access details: https://github.com/JBEI/ART/#license
- Source repository: https://github.com/JBEI/AutomatedRecommendationTool

Because ART requires a license to run, think about who can reach the machine and
directory where you place the code before you clone it.

## 2. Choose how to run ART

There are two routes. Pick based on how you want ART to sit next to SynThera.

| Route | What you get | When to use it |
|---|---|---|
| A. ART's own Docker image | ART running in its own `jbei/art` container with its own Jupyter server and tutorials | You want the environment ART's maintainers build and support |
| B. ART inside the SynThera environment | ART importable on the SynThera kernel, so `Step2_BuildModel` runs alongside the rest of the pipeline | You want the whole SynThera DBTL loop to run in one place |

ART's maintainers support Route A and do not support direct installs into another
Python environment. Route B is still valid (the SynThera authors use it) and
mirrors what ART's own image does internally, but you own the setup if a
dependency shifts.

## Route A: run ART's own Docker image (maintainer-supported)

From the ART install guide, with Docker installed and running:

```bash
# clone into a directory named 'art' (the container takes this name)
git clone https://github.com/JBEI/AutomatedRecommendationTool.git art
cd art

# build the licensed image locally (it is not on Docker Hub)
docker build --pull -t jbei/art .

# generate compose.yml, environment.env, and run_jupyter.sh in the directory
docker run -it --rm -v "$(pwd)":/app jbei/art create_config_files

# launch ART's Jupyter server (local only) and open it in the browser
./run_jupyter.sh
```

Notebooks you want ART to run must live in this install directory. Verify the
install by running `notebooks/Limonene_Example.ipynb`. This route keeps ART in
its own container, separate from the SynThera image.

## Route B: install ART into the SynThera environment (integrated)

This makes `import art` work on the SynThera kernel, so `Step2_BuildModel` runs
in the same container as the other steps. You still need the license and the ART
source from step 1. The SynThera build exposes two mechanisms for this.

**At build time**, pass the ART source as a build argument (run from `build/`):

```bash
docker build \
  --build-arg ART_SOURCE="git+https://github.com/JBEI/AutomatedRecommendationTool.git" \
  -t synthera:with-art .
```

Or uncomment the `args` block in `docker-compose.yml` and set the same value,
then `docker compose up --build`.

**Into a running container**, use the helper script, then commit so the change
survives a restart:

```bash
./installArt.sh synthera "git+https://github.com/JBEI/AutomatedRecommendationTool.git"
docker commit synthera synthera:with-art
```

If a plain `pip install` of the source does not resolve, that is expected given
ART is Docker-first. Fall back to a local checkout: clone the repo, follow ART's
own `Dockerfile` as the reference for its dependencies, build a wheel or install
it editable inside the container with `docker exec`, then `docker commit`. If ART
needs its configuration files or thread settings (its `create_config_files` step
and the `OPENBLAS_NUM_THREADS` / `OMP_NUM_THREADS` values it sets to `1`),
replicate those in the SynThera container the same way before committing.

## 3. Add ART to the workflow

Once ART is installed by Route B:

- Select the **SynThera** kernel and run `Step2_BuildModel.ipynb`; `import art`
  now succeeds.
- The other steps (1, 3, 4, and Step 5 on `SynThera_ExptDesigns`) do not use ART
  and are unaffected.

Two settings matter for reproducible ART runs, since ART uses all available
cores by default: set a fixed `seed`, and cap `max_mcmc_cores`. On a shared HPC
node, also set `OPENBLAS_NUM_THREADS` and `OMP_NUM_THREADS` to control low-level
threading, as ART's own image does.
