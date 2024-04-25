#! /usr/bin/env bash

mkdir -p ${CONDA_PREFIX}/bin

# compilation
g++ src/harpy/globalscripts/extractReads.cpp -O3 -o ${CONDA_PREFIX}/bin/extractReads

## associated scripts
chmod +x src/harpy/globalscripts/*
cp src/harpy/globalscripts/* ${CONDA_PREFIX}/bin/

