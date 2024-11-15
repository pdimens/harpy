#! /usr/bin/env bash

## Use the first positional argument to set a name, usually `harpy` or `harpytest`

if command -v conda &> /dev/null
then
    conda env create -n $1 -f resources/harpy.yaml
else
    conda env create -n $1 -f resources/harpy.yaml
fi
