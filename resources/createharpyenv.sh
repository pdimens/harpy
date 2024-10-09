#! /usr/bin/env bash

## Use the first positional argument to set a name, usually `harpy` or `harpytest`

mamba create -n $1 -f harpy.yaml