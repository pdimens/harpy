#! /usr/bin/env bash

mkdir -p ${CONDA_PREFIX}/bin

# rules
cp rules/*.smk ${CONDA_PREFIX}/bin/

# associated scripts
chmod +x utilities/*
cp utilities/* ${CONDA_PREFIX}/bin/

# reports
cp reports/*.Rmd ${CONDA_PREFIX}/bin/
