#! /usr/bin/env bash

mkdir -p ${CONDA_PREFIX}/bin

# rules
cp workflow/rules/*.smk ${CONDA_PREFIX}/bin/

# compilation
g++ workflow/scripts/extractReads.cpp -O3 -o workflow/scripts/extractReads && rm workflow/scripts/extractReads.cpp

# associated scripts
chmod +x workflow/scripts/*
cp workflow/scripts/* ${CONDA_PREFIX}/bin/

# reports
cp workflow/report/*.Rmd ${CONDA_PREFIX}/bin/