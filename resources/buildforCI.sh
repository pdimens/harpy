#! /usr/bin/env bash

mkdir -p ${CONDA_PREFIX}/bin

# compilation
g++ src/harpy/bin/extractReads.cpp -O3 -o ${CONDA_PREFIX}/bin/extractReads

## associated scripts
chmod +x src/harpy/bin/*
cp src/harpy/bin/* ${CONDA_PREFIX}/bin/

