#! /usr/bin/env bash
set -euo pipefail 

if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: active conda environment not detected."
    echo "To use this installation script, you need to already be in an active conda environment."
    exit 1
fi

mkdir -p ${CONDA_PREFIX}/bin

# install harpy proper
pip install --no-deps --disable-pip-version-check -e . && rm -rf build

# shell completion scripts (bash, zsh, fish), sourced on activation
bash resources/generate_completion.sh "${CONDA_PREFIX}" python
mkdir -p "${CONDA_PREFIX}/etc/conda/activate.d"
cp resources/shell_completion.sh "${CONDA_PREFIX}/etc/conda/activate.d/harpy-completion.sh"
cp resources/shell_completion.fish "${CONDA_PREFIX}/etc/conda/activate.d/harpy-completion.fish"

go build -C harpy/utils/demoxolotl -o $CONDA_PREFIX/bin/demoxolotl -ldflags='-s -w'
go build -C harpy/utils/standardize -o $CONDA_PREFIX/bin/djinn-standardize -ldflags='-s -w' standardize.go 
