#! /usr/bin/env bash

if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: active conda environment not detected."
    echo "To use this installation script, you need to already be in an active conda environment."
    exit 1
fi

mkdir -p ${CONDA_PREFIX}/bin

# install harpy proper
pip install . --no-deps

# rules
cp -f workflow/rules/*.smk ${CONDA_PREFIX}/bin/

# associated scripts
chmod +x workflow/scripts/* 
cp -f workflow/scripts/* ${CONDA_PREFIX}/bin/

# reports
cp -f workflow/report/*.Rmd ${CONDA_PREFIX}/bin/

# completion
#cp misc/harpy_completion.sh ${CONDA_PREFIX}/etc/conda/activate.d/