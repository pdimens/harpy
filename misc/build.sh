#! /usr/bin/env bash

# install LepWrap into conda PATH
mkdir -p $CONDA_PREFIX/bin

# build and install ema
cd ema
make
chmod +x ema
cp ema $CONDA_PREFIX/bin/ema-h
cd ..

# Harpy executable
chmod +x harpy
cp harpy $CONDA_PREFIX/bin/

# rules
cp rules/*.smk $CONDA_PREFIX/bin/

# associated scripts
chmod +x utilities/*.{py,r,pl}
cp utilities/*.{py,r,pl} $CONDA_PREFIX/bin/

# reports
chmod +x reports/*.Rmd
cp reports/*.Rmd $CONDA_PREFIX/bin/