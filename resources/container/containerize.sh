#! /usr/bin/env bash

snakemake -s containerize.smk --containerize > Dockerfile

docker build -t pdimens/harpy:1.0 .