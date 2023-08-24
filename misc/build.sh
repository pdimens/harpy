#! /usr/bin/env bash

mkdir -p ${PREFIX}/bin
#cp misc/ema-h ${PREFIX}/bin

# Harpy executable
#cp harpy ${PREFIX}/bin/
#pip install -e .

# rules
cp rules/*.smk ${PREFIX}/bin/

# associated scripts
cp utilities/* ${PREFIX}/bin/

# reports
cp reports/*.Rmd ${PREFIX}/bin/