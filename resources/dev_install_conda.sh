#! /usr/bin/env bash

conda env create --prefix harpy/.conda/harpy $1 --file harpy/resources/harpy.yaml

conda activate harpy/.conda/harpy

mkdir -p ${CONDA_PREFIX}/bin

cd harpy

# compilation
g++ harpy/bin/extractReads.cpp -O3 -o ${CONDA_PREFIX}/bin/extractReads

# install harpy proper
pip install --no-deps --disable-pip-version-check -e . && \
    rm -rf build

# associated scripts
chmod +x harpy/bin/* 
cp -f harpy/bin/* ${CONDA_PREFIX}/bin/