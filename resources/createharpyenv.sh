#! /usr/bin/env bash

## Use the first positional argument to set a name, usually `harpy` or `harpytest`

mamba create -n $1 -c bioconda -c conda-forge \
    bcftools=1.19 \
    python \
    pandas \
    rich-click \
    snakemake-minimal>7 \
    samtools \
    seqtk