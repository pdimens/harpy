#! /usr/bin/env bash

mkdir -p ${PREFIX}/bin

# compilation
g++ src/harpy/scripts/extractReads.cpp -O3 -o ${PREFIX}/bin/extractReads

# install harpy proper
${PYTHON} -m pip install . --no-deps -vvv

# associated scripts
chmod +x globalscripts/* 
cp globalscripts/* ${PREFIX}/bin/
