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
cp -f src/harpy/rules/*.smk ${CONDA_PREFIX}/bin/

# compilation
g++ src/harpy/scripts/extractReads.cpp -O3 -o src/harpy/scripts/extractReads 
#&& rm src/harpy/scripts/extractReads.cpp

#g++ -O3 -o src/harpy/scripts/demuxGen1 src/harpy/scripts/demult_fastq.cpp -lgzstream -I/usr/local/include/gzstream -I/usr/include -std=gnu++11 -lz -Wall

# associated scripts
chmod +x src/harpy/scripts/* 
cp -f src/harpy/scripts/* ${CONDA_PREFIX}/bin/

# reports
cp -f src/harpy/reports/*.Rmd ${CONDA_PREFIX}/bin/