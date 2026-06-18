#!/usr/bin/env bash
# bootstrap.sh — one command to make a fresh machine MuAgene-ready.
#
# Stands up the single integrated `muagene` conda env (science stack + both agent CLIs) and
# hands off to `Execution-MuAgent init-machine`, which validates the env, installs the agent
# packages, writes ~/.muagene/machine.config, and records the env fingerprint. After it prints
# "Machine ready.", `conda activate muagene` is all that's needed to drive the agents.
#
# There is deliberately NO separate bootstrap env: one env serves the science workflow, the
# interactive `Processing-MuAgent`/`Execution-MuAgent` CLIs, and the `python -m execution_muagent`
# daemon that `Processing-MuAgent submit` spawns.
#
# Usage:
#   bash Execution-MuAgent/scripts/bootstrap.sh \
#        [--processing-repo DIR] [--device cpu|gpu|both] [--conda-env NAME] \
#        [extra init-machine flags, e.g. --gpu-image-uri docker://reg/muagene-gpu:TAG]
#
# Defaults: --processing-repo = the sibling ../Processing-MuAgent of this repo; --device cpu.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXEC_REPO="$(cd "$SCRIPT_DIR/.." && pwd)"                       # Execution-MuAgent checkout
PROCESSING_REPO="$(cd "$EXEC_REPO/.." && pwd)/Processing-MuAgent"
DEVICE="cpu"
ENV_NAME="muagene"
PASSTHRU=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --processing-repo) PROCESSING_REPO="$2"; shift 2;;
    --device)          DEVICE="$2"; shift 2;;
    # init-machine owns the env name; mirror it here so we create the same env, and pass through.
    --conda-env)       ENV_NAME="$2"; PASSTHRU+=("$1" "$2"); shift 2;;
    *)                 PASSTHRU+=("$1"); shift;;
  esac
done

# 1. Resolve a conda env manager (prefer micromamba > mamba > conda).
MANAGER=""
for m in micromamba mamba conda; do
  if command -v "$m" >/dev/null 2>&1; then MANAGER="$m"; break; fi
done
if [[ -z "$MANAGER" ]]; then
  echo "ERROR: no conda env manager (micromamba/mamba/conda) on PATH. Install miniforge first." >&2
  exit 1
fi

# 2. The CPU env installs from the committed, solve-free lock — fail loud if it's missing.
LOCK="$PROCESSING_REPO/workflow/envs/processing.linux-64.lock"
if [[ ! -f "$LOCK" ]]; then
  echo "ERROR: CPU lock not found at $LOCK — is --processing-repo a Processing-MuAgent checkout?" >&2
  exit 1
fi

echo "==> manager=$MANAGER  processing_repo=$PROCESSING_REPO  device=$DEVICE  env=$ENV_NAME"

# 3. Create the integrated env from the lock if absent. Creating it HERE (not inside
#    init-machine) means init-machine never mutates the interpreter it runs on — it only
#    validates, records the fingerprint, and installs the editable agent packages.
if "$MANAGER" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "==> env '$ENV_NAME' already present; init-machine will reconcile it"
else
  echo "==> creating '$ENV_NAME' from $(basename "$LOCK") (solve-free)"
  "$MANAGER" create -y -n "$ENV_NAME" --file "$LOCK"
fi

# 4. Make the Execution CLI runnable inside the env so we can invoke init-machine from it.
#    init-machine re-affirms this (and adds Processing) with the same --no-deps editable install.
echo "==> installing Execution-MuAgent CLI into '$ENV_NAME'"
"$MANAGER" run -n "$ENV_NAME" pip install --no-deps -e "$EXEC_REPO"

# 5. Hand off to the provisioning engine: validate, install both agent packages, write the
#    machine profile, record the fingerprint, and (for --device gpu/both) pull the GPU image.
echo "==> running init-machine"
exec "$MANAGER" run -n "$ENV_NAME" Execution-MuAgent init-machine \
  --processing-repo "$PROCESSING_REPO" --device "$DEVICE" ${PASSTHRU[@]+"${PASSTHRU[@]}"}
