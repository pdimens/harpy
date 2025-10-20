#! /usr/bin/env python

import shutil
import subprocess
import os

environ = {
    "align" : [
        "bwa-mem2",
        "bwa",
        "samtools==1.22",
        "seqtk",
        "strobealign",
        "tabix"
    ],
    "assembly" : [
        "arcs",
        "bwa",
        "cloudspades",
        "links",
        "quast",
        "busco",
        "samtools",
        "tigmint"
    ],
    "deconvolution" : [
        "quickdeconvolution"
    ],

    "demultiplex": [
        "dmox>=0.2"
    ],
    "metassembly": [
        "athena_meta==1.2"
    ],
    "phase" : [
        "hapcut2",
        "whatshap"
    ],
    "qc" : [
        "click==8.2.1",
        "falco==1.2.5",
        "fastp",
        "multiqc==1.30",
        "pysam==0.23"
    ],
    "report" : [
        "quarto",
        "r-dt",
        "r-dplyr",
        "r-highcharter",
        "r-magrittr",
        "r-plotly",
        "r-scales",
        "r-tidyr",
        "r-viridislite", 
        "r-xml2",
        "r-biocircos"
    ],
    "simulations" : [
        "simug>=1.0.1"
    ],
    "stitch" : [
        "r-stitch>=1.8.4"
    ],
    "variants" : [
        "bcftools==1.22",
        "freebayes==1.3.9",
        "leviathan",
        "naibr-plus",
        "setuptools"
    ]
}

dockerfile_text = """
FROM ghcr.io/prefix-dev/pixi:0.56.0 AS build

# copy source code, pixi.toml and pixi.lock to the container
WORKDIR /app
COPY . .

# use `--locked` to ensure the lockfile is up to date with pixi.toml
RUN pixi install --locked && rm -rf ~/.cache/rattler

# create the shell-hook bash script to activate the environment
RUN echo "#!/bin/bash" > /app/entrypoint.sh && \\
    pixi shell-hook -s bash >> /app/entrypoint.sh && \\
    echo 'exec "$@"' >> /app/entrypoint.sh && \\
    chmod +x /app/entrypoint.sh

FROM ubuntu:24.04 AS production
WORKDIR /app
COPY --from=build --chmod=0755 /app/entrypoint.sh /app/entrypoint.sh
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default

ENTRYPOINT ["/app/entrypoint.sh"]
"""

def create_pixi_dockerfiles():
    '''
    Using the defined environments, create a series of folders where each has a dockerfile
    and pixi.toml file to create one of the environments.
    '''
    shutil.rmtree("container", ignore_errors=True)
    for env,deps in environ.items():
        os.makedirs(f"container/{env}", exist_ok=True)
        with open(f"container/{env}/Dockerfile", "w") as dockerfile:
            dockerfile.write(dockerfile_text)
        if env == "report":
            subprocess.run(f"pixi init container/{env} -c conda-forge -c r".split())
        else:
            subprocess.run(f"pixi init container/{env} -c conda-forge -c bioconda".split())

        subprocess.run(
            ["pixi", "add", "--no-progress", "--manifest-path", f"container/{env}/pixi.toml", *deps]
        )
        shutil.rmtree("container/.pixi", ignore_errors=True)
