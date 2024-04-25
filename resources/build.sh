#! /usr/bin/env bash

mkdir -p ${PREFIX}/bin

# compilation
g++ src/harpy/globalscripts/extractReads.cpp -O3 -o ${PREFIX}/bin/extractReads

# install harpy proper
${PYTHON} -m pip install . --no-deps -vvv

# associated scripts
chmod +x src/harpy/globalscripts/* 
cp src/harpy/globalscripts/* ${PREFIX}/bin/
