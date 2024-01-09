#! /usr/bin/env bash

docker build -t qc.fat --target qc -f $PWD/resources/container/qc/Dockerfile $PWD

# minify the docker image

sudo slim build --http-probe=false --include-shell --tag qc qcimage