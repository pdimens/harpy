#! /usr/bin/env bash
# Writes harpy's tab-completion scripts for bash, zsh, and fish to <prefix>/share/harpy/
# (complete.bash, complete.zsh, complete.fish), where the activation hooks in
# shell_completion.sh and shell_completion.fish expect to find them.
#
# usage: generate_completion.sh [PREFIX] [PYTHON]
#   PREFIX  environment to write to, must already have harpy installed (default: $CONDA_PREFIX)
#   PYTHON  python interpreter of that environment (default: python)
set -euo pipefail

PREFIX_DIR="${1:-${CONDA_PREFIX:?Error: no PREFIX given and no active conda/pixi environment detected.}}"
PYTHON_BIN="${2:-python}"

mkdir -p "${PREFIX_DIR}/share/harpy"
for shell in bash zsh fish; do
    "${PYTHON_BIN}" -m harpy completion "${shell}" > "${PREFIX_DIR}/share/harpy/complete.${shell}"
done
