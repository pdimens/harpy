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
    cd harpy/utils/stagger && go build -ldflags="-s -w" -o ../gih-stagger stagger.go
    cd ../convert && go build -ldflags="-s -w" -o ../gih-convert convert.go
    cd ../standardize && go build -ldflags="-s -w" -o ../djinn-standardize standardize.go
    cd ..
    chmod +x gih-stagger gih-convert djinn-standardize
    mv gih-stagger gih-convert djinn-standardize ${CONDA_PREFIX}/bin/
    cd ../..
}

# associated scripts
#chmod +x harpy/bin/* 
#cp -f harpy/bin/* ${CONDA_PREFIX}/bin/