#! /usr/bin/env bash

## Use the first positional argument to set a name, usually `harpy` or `harpytest`

conda env create -n $1 -f resources/harpy.yaml