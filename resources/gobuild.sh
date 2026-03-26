#! /usr/bin/env bash

{
    cd harpy/utils/stagger
    go mod tidy && go build -ldflags="-s -w" -o ../gih-stagger stagger.go
    cd ../convert
    go mod tidy && go build -ldflags="-s -w" -o ../gih-convert convert.go
    cd ../standardize
    go mod tidy && go build -ldflags="-s -w" -o ../djinn-standardize standardize.go
    cd .. && chmod +x gih-stagger gih-convert djinn-standardize
    mv gih-stagger gih-convert djinn-standardize ${CONDA_PREFIX}/bin/
}
