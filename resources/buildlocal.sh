#! /usr/bin/env bash

if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: active conda environment not detected."
    echo "To use this installation script, you need to already be in an active conda environment."
    exit 1
fi

mkdir -p ${CONDA_PREFIX}/bin

# install harpy proper
pip install --no-deps --disable-pip-version-check -e . && rm -rf build

{
    cd harpy/utils/preproc
    go build -ldflags="-s -w" -o gih-stagger stagger.go
    go build -ldflags="-s -w" -o gih-convert convert/convert.go
    chmod +x gih-stagger gih-convert
    mv gih-stagger gih-convert ${CONDA_PREFIX}/bin/
    cd ../../..
}

# associated scripts
#chmod +x harpy/bin/* 
#cp -f harpy/bin/* ${CONDA_PREFIX}/bin/