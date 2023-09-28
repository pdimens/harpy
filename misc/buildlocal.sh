#! /usr/bin/env bash

if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: active conda environment not detected."
    echo "To use this installation script, you need to already be in an active conda environment."
    exit 1
fi

mkdir -p ${CONDA_PREFIX}/bin

# install harpy proper
pip install . --ignore-installed --no-deps

# rules
cp -f rules/*.smk ${CONDA_PREFIX}/bin/

# associated scripts
chmod +x utilities/* 
cp -f utilities/* ${CONDA_PREFIX}/bin/

# reports
cp -f reports/*.Rmd ${CONDA_PREFIX}/bin/
