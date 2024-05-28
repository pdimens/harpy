#! /usr/bin/env bash

snakemake -s resources/container/containerize.smk --directory . --containerize > Dockerfile

docker build -t pdimens/harpy .

mv Dockerfile resources/container/