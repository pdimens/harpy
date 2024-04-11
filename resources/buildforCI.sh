#! /usr/bin/env bash

mkdir -p ${CONDA_PREFIX}/bin

# compilation
g++ src/harpy/scripts/extractReads.cpp -O3 -o ${CONDA_PREFIX}/bin/extractReads

## associated scripts
chmod +x globalscripts/*
cp globalscripts/* ${CONDA_PREFIX}/bin/

