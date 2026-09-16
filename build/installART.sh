#!/usr/bin/env bash
# Install ART (Automated Recommendation Tool) into the SynThera environment of a
# running SynThera container, so Step2_BuildModel can import it.
#
# ART is free for non-commercial academic, non-profit, and U.S. government use
# but licensed separately, and is intentionally not part of the image. Obtain a
# license and the source first (see ART.md and https://github.com/JBEI/ART/#license).
# ART's maintainers support their own Docker image rather than direct installs,
# so if the pip install below does not resolve, follow ART.md's fallback.
#
# The art source you pass must be reachable from inside the container: a git URL,
# a pip spec, or a wheel/source path that exists in the container filesystem.

set -euo pipefail

containerName="${1:-synthera}"
artSource="${2:-}"

if [ -z "$artSource" ]; then
  echo "Usage: ./installArt.sh <container_name> <art_source>"
  echo
  echo "  container_name  name of the running container (default: synthera)"
  echo "  art_source      git URL, pip spec, or in-container wheel/source path"
  echo
  echo "Example (after obtaining a license):"
  echo "  ./installArt.sh synthera \\"
  echo "    'git+https://github.com/JBEI/AutomatedRecommendationTool.git'"
  exit 1
fi

echo "Installing ART into SynThera in container '${containerName}'."
docker exec "${containerName}" conda run -n SynThera pip install "${artSource}"

echo
echo "ART installed. Commit the container so the change survives a restart:"
echo "  docker commit ${containerName} synthera:with-art"
echo
echo "Then select the SynThera kernel and run Step2_BuildModel.ipynb."
