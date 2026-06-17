#!/usr/bin/env bash
# build_and_push_gpu_image.sh — MAINTAINER-ONLY, run OUT OF BAND (not part of machine
# setup). Builds the MuAgene GPU container from workflow/envs/muagene-gpu.def and pushes
# it to a registry as a pinned image. Target machines then PULL that image
# (Execution-MuAgent provision-env/init-machine --gpu-image-uri); no machine builds the
# container locally — that is the whole point of pull-only (no --fakeroot/subuid needed
# on compute hosts).
#
# Run this once per change to muagene-gpu.def, on a host that CAN build a container:
#   - a box where your user is in /etc/subuid (unprivileged --fakeroot build), OR
#   - `singularity build --remote` with a Sylabs token, OR
#   - a docker/buildah host.
# Then BUMP the published tag — republishing a new tag is what flips every machine's
# GPU env to `stale` (the fingerprint is the pinned image reference).
#
# Usage:
#   bash scripts/build_and_push_gpu_image.sh <image_uri> [singularity_module]
# Example:
#   bash scripts/build_and_push_gpu_image.sh oras://registry.example.org/muagene-gpu:25.04 \
#        singularityce/3.11.3
set -euo pipefail

IMAGE_URI="${1:-}"
SINGULARITY_MODULE="${2:-}"
if [ -z "$IMAGE_URI" ]; then
    echo "usage: $0 <image_uri> [singularity_module]" >&2
    echo "  <image_uri> e.g. oras://<registry>/muagene-gpu:<tag> or docker://<registry>/muagene-gpu:<tag>" >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DEF="$REPO_ROOT/workflow/envs/muagene-gpu.def"
[ -f "$DEF" ] || { echo "missing recipe: $DEF" >&2; exit 1; }

if [ -n "$SINGULARITY_MODULE" ]; then
    # shellcheck disable=SC1091
    module load "$SINGULARITY_MODULE" 2>/dev/null || true
fi

RUNTIME="$(command -v apptainer || command -v singularity || true)"
[ -n "$RUNTIME" ] || { echo "no apptainer/singularity on PATH" >&2; exit 1; }

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
SIF="$WORK/muagene-gpu.sif"

echo ">> building $SIF from $DEF"
# --fakeroot for the %post pip layer; needs userns/subuid on THIS build host (this is
# the only place a build happens — never on a compute/target machine).
"$RUNTIME" build --fakeroot "$SIF" "$DEF"

echo ">> pushing $SIF -> $IMAGE_URI"
"$RUNTIME" push "$SIF" "$IMAGE_URI"

echo ">> published $IMAGE_URI"
echo "   Point machines at it:  Execution-MuAgent init-machine --device gpu --gpu-image-uri $IMAGE_URI"
echo "   (or per run: Processing-MuAgent configure-execution ... --gpu-image-uri $IMAGE_URI)"
