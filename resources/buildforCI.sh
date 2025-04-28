#! /usr/bin/env bash

mkdir -p ${CONDA_PREFIX}/bin

## associated scripts
chmod +x harpy/bin/*
cp harpy/bin/* ${CONDA_PREFIX}/bin/

