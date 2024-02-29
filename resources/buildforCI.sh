#! /usr/bin/env bash

mkdir -p ${CONDA_PREFIX}/bin

# rules
cp workflow/rules/*.smk ${CONDA_PREFIX}/bin/

# associated scripts
chmod +x workflow/scripts/*
cp workflow/scripts/* ${CONDA_PREFIX}/bin/

# reports
cp workflow/report/*.Rmd ${CONDA_PREFIX}/bin/