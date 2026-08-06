#! /usr/bin/env bash

python -m ipykernel install --prefix "${CONDA_PREFIX:?CONDA_PREFIX is required}" --name ipython-harpy \
    --display-name "Python (harpy)"