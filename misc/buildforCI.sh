#! /usr/bin/env bash

mkdir -p ${CONDA_PREFIX}/bin
#cp -n misc/ema-h ${CONDA_PREFIX}/bin

# Harpy executable
#cp harpy ${CONDA_PREFIX}/bin/

# rules
cp rules/*.smk ${CONDA_PREFIX}/bin/

# associated scripts
cp utilities/* ${CONDA_PREFIX}/bin/

# reports
cp reports/*.Rmd ${CONDA_PREFIX}/bin/
