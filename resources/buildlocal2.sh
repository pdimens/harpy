#! /usr/bin/env bash

if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: active conda environment not detected."
    echo "To use this installation script, you need to already be in an active conda environment."
    exit 1
fi

mkdir -p ${CONDA_PREFIX}/bin

# compilation
g++ src/harpy/scripts/extractReads.cpp -O3 -o ${CONDA_PREFIX}/bin/extractReads

# install harpy proper
pip install . --no-deps

## rules
#cp -f workflow/rules/*.smk ${CONDA_PREFIX}/bin/

#&& rm workflow/scripts/extractReads.cpp

#g++ -O3 -o workflow/scripts/demuxGen1 workflow/scripts/demult_fastq.cpp -lgzstream -I/usr/local/include/gzstream -I/usr/include -std=gnu++11 -lz -Wall

# associated scripts
#chmod +x workflow/scripts/* 
#cp -f workflow/scripts/* ${CONDA_PREFIX}/bin/

# reports
#cp -f workflow/report/*.Rmd ${CONDA_PREFIX}/bin/