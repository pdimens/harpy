#! /usr/bin/env bash
# Writes harpy's tab-completion scripts for bash, zsh, and fish into <prefix>/share/ at the
# standard locations each shell searches for completions:
#   share/bash-completion/completions/harpy
#   share/zsh/site-functions/_harpy
#   share/fish/vendor_completions.d/harpy.fish
# The activation hooks (shell_completion.sh / .fish) make the shells look there.
#
# usage: generate_completion.sh [PREFIX] [PYTHON]
#   PREFIX  environment to write to, must already have harpy installed (default: $CONDA_PREFIX)
#   PYTHON  python interpreter of that environment (default: python)
set -euo pipefail

PREFIX_DIR="${1:-${CONDA_PREFIX:?Error: no PREFIX given and no active conda/pixi environment detected.}}"
PYTHON_BIN="${2:-python}"

mkdir -p "${PREFIX_DIR}/share/bash-completion/completions" "${PREFIX_DIR}/share/zsh/site-functions" "${PREFIX_DIR}/share/fish/vendor_completions.d"
"${PYTHON_BIN}" -m harpy completion bash > "${PREFIX_DIR}/share/bash-completion/completions/harpy"
"${PYTHON_BIN}" -m harpy completion zsh  > "${PREFIX_DIR}/share/zsh/site-functions/_harpy"
"${PYTHON_BIN}" -m harpy completion fish > "${PREFIX_DIR}/share/fish/vendor_completions.d/harpy.fish"
