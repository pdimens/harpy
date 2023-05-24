#! /usr/bin/env bash

mkdir -p ${CONDA_PREFIX}/bin
cp misc/ema-h ${CONDA_PREFIX}/bin

# Harpy executable
cp harpy ${CONDA_PREFIX}/bin/

# rules
cp rules/*.smk ${CONDA_PREFIX}/bin/

# associated scripts
cp utilities/*.{py,R,pl} ${CONDA_PREFIX}/bin/

# reports
cp reports/*.Rmd ${CONDA_PREFIX}/bin/


### deprecated? ###
#export C_INCLUDE_PATH=${PREFIX}/include
#export CPLUS_INCLUDE_PATH=${PREFIX}/include
#export CPP_INCLUDE_PATH=${PREFIX}/include
#export CXX_INCLUDE_PATH=${PREFIX}/include
#export LIBRARY_PATH=${PREFIX}/lib

# create conda PATH, if not already existing

# build and install ema
#rm -rf ema || true
#git clone --recursive https://github.com/EdHarry/ema.git
#cd ema && git checkout haplotag
#git submodule update --remote
#git apply ../harpy/misc/makefile.patch
#git apply ../harpy/misc/remove_native.patch
#make CC="$CC -fcommon" CXX="$CXX -fcommon" LDFLAGS="$LDFLAGS"
#chmod +x ema
#cp ema $CONDA_PREFIX/bin/ema-h
#cd ..