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
    cd harpy/utils
    go build -C stagger -o ../gih-stagger -ldflags='-s -w' stagger.go
    go build -C convert -o ../gih-convert -ldflags='-s -w' convert.go
    go build -C standardize -o ../djinn-standardize -ldflags='-s -w' standardize.go 
    chmod +x gih-stagger gih-convert djinn-standardize
    mv gih-stagger gih-convert djinn-standardize "${CONDA_PREFIX}/bin/"
    cd ../..
}
