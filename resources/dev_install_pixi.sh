#! /usr/bin/env bash

ENV_PREFIX="harpy/.pixi/envs/default/bin/"

if ! pixi init -i harpy/resources/harpy.yaml harpy; then
    echo "Error: Failed to initialize pixi environment"
    exit 1
fi

if ! echo -e "\n[pypi-dependencies]\nharpy = { path = \".\", editable = true}" >> harpy/pixi.toml; then
    echo "Error: Failed to update pixi.toml"
    exit 1
fi

if ! cd harpy; then
    echo "Error: Failed to change directory to harpy"
    exit 1
fi

pixi install
