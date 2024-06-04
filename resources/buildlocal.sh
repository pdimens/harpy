#! /usr/bin/env bash

if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: active conda environment not detected."
    echo "To use this installation script, you need to already be in an active conda environment."
    exit 1
fi

mkdir -p ${CONDA_PREFIX}/bin

# compilation
g++ src/harpy/bin/extractReads.cpp -O3 -o ${CONDA_PREFIX}/bin/extractReads

# install harpy proper
pip install . --no-deps

# associated scripts
chmod +x src/harpy/bin/* 
cp -f src/harpy/bin/* ${CONDA_PREFIX}/bin/