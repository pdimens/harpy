#! /usr/bin/env bash

mkdir -p ${PREFIX}/bin

# install harpy proper
${PYTHON} -m pip install . --no-deps -vvv

# rules
cp workflow/rules/*.smk ${PREFIX}/bin/

# compilation
g++ workflow/scripts/extractReads.cpp -O3 -o workflow/scripts/extractReads && rm workflow/scripts/extractReads.cpp

# associated scripts
chmod +x workflow/scripts/* 
cp workflow/scripts/* ${PREFIX}/bin/

# reports
cp workflow/report/*.Rmd ${PREFIX}/bin/