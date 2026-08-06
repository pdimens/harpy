#! /usr/bin/env bash

python -m ipykernel install --prefix "$CONDA_PREFIX" --name ipython-harpy \
    --display-name "Python (harpy)" >/dev/null 2>&1