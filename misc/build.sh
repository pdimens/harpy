#! /usr/bin/env bash

mkdir -p ${PREFIX}/bin
#cp misc/ema-h ${PREFIX}/bin

# install harpy proper
${PREFIX}/bin/python -m pip install . --no-deps -vv

# rules
cp rules/*.smk ${PREFIX}/bin/

# associated scripts
chmod +x utilities/* 
cp utilities/* ${PREFIX}/bin/

# reports
cp reports/*.Rmd ${PREFIX}/bin/

# completion
cp misc/harpy_completion.sh ${PREFIX}/etc/conda/activate.d/