#! /usr/bin/env bash

ENV_PREFIX="harpy/.pixi/envs/default/bin/"

pixi init -i harpy/resources/harpy.yaml harpy
echo -e "\n[pypi-dependencies]\nharpy = { path = \".\", editable = true}" >> harpy/pixi.toml

cd harpy && pixi shell

mkdir -p ${ENV_PREFIX}/bin

# compilation
g++ harpy/bin/extractReads.cpp -O3 -o ${ENV_PREFIX}/bin/extractReads

# associated scripts
chmod +x harpy/bin/* 
cp -f harpy/bin/* ${ENV_PREFIX}/bin/
