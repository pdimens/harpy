#! /usr/bin/env bash

if ! conda env create --prefix harpy/.conda/harpy $1 --file harpy/resources/harpy.yaml; then
    echo "Error: Failed to create conda environment"
    exit 1
fi

if ! conda activate harpy/.conda/harpy; then
cd harpy || { echo "Error: Failed to change directory to harpy"; exit 1; }

fi
mkdir -p ${CONDA_PREFIX}/bin

# install harpy proper
if ! pip install --no-deps --disable-pip-version-check -e .; then
    echo "Error: Failed to install harpy package"
    exit 1
fi

# install harpy proper
pip install --no-deps --disable-pip-version-check -e . && \
    rm -rf build

# associated scripts
chmod +x harpy/bin/* 
cp -f harpy/bin/* ${CONDA_PREFIX}/bin/