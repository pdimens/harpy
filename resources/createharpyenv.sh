#! /usr/bin/env bash

## Use the first positional argument to set a name, usually `harpy` or `harpytest`

if ! command -v conda &> /dev/null; then
    echo "Error: conda is not installed or not in PATH" >&2
    exit 1
fi

conda env create -n $1 -f resources/harpy.yaml