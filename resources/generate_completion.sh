#! /usr/bin/env bash
# Writes the tab-completion scripts for harpy, harpy-utils, and hv for bash, zsh, and fish into <prefix>/share/ at the
# standard locations each shell searches for completions:
#   share/bash-completion/completions/<program>
#   share/zsh/site-functions/_<program>
#   share/fish/vendor_completions.d/<program>.fish
# The activation hooks (shell_completion.sh / .fish) make the shells look there.
#
# usage: generate_completion.sh [PREFIX] [PYTHON]
#   PREFIX  environment to write to, must already have harpy installed (default: $CONDA_PREFIX)
#   PYTHON  python interpreter of that environment (default: python)
set -euo pipefail

PREFIX_DIR="${1:-${CONDA_PREFIX:?Error: no PREFIX given and no active conda/pixi environment detected.}}"
PYTHON_BIN="${2:-python}"

mkdir -p "${PREFIX_DIR}/share/bash-completion/completions" "${PREFIX_DIR}/share/zsh/site-functions" "${PREFIX_DIR}/share/fish/vendor_completions.d"
for program in harpy harpy-utils hv; do
    "${PYTHON_BIN}" -m harpy completion bash "${program}" > "${PREFIX_DIR}/share/bash-completion/completions/${program}"
    "${PYTHON_BIN}" -m harpy completion zsh  "${program}" > "${PREFIX_DIR}/share/zsh/site-functions/_${program}"
    "${PYTHON_BIN}" -m harpy completion fish "${program}" > "${PREFIX_DIR}/share/fish/vendor_completions.d/${program}.fish"
done
