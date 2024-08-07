#! /usr/bin/env bash

mkdir -p ${CONDA_PREFIX}/bin

# compilation
g++ harpy/bin/extractReads.cpp -O3 -o ${CONDA_PREFIX}/bin/extractReads

## associated scripts
chmod +x harpy/bin/*
cp harpy/bin/* ${CONDA_PREFIX}/bin/

