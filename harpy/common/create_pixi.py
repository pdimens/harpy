#! /usr/bin/env python

import glob
import shutil
import subprocess
import os
import sys

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

def create_pixi_dockerfiles():
    '''
    Using the defined environments, create a series of folders where each has a dockerfile to create one of the environments.
    '''
    rm_cache = "&& rm -rf home/.cache/rattler".split()
    for env,deps in environ.items():
        os.makedirs(f"container/{env}", exist_ok=True)
        with open(f"container/{env}/Dockerfile", "w") as dockerfile:
            dockerfile.write("FROM ghcr.io/prefix-dev/pixi:bookworm-slim\n\nRUN ")
            if env == "report":
                channels = ["conda-forge", "r"]
            else:
                channels = ["conda-forge", "bioconda"]
            _chan = " ".join([f"--channel {i}" for i in channels]).split()
            if env == "deconvolution":
                dockerfile.write(
                    " ".join(["pixi", "global", "install", *_chan, "--environment", env, "--expose", "QuickDeconvolution", *deps, *rm_cache])
                )
            else:
                dockerfile.write(
                    " ".join(["pixi", "global", "install", *_chan, "--environment", env, *deps, *rm_cache])
                )


def reset_pixi_global():
    # remove any existing global packages
    pixidir = os.environ['HOME'] + "/.pixi"
    for f in glob.glob(pixidir + "/bin/*"):
        if os.path.isdir(f):
            shutil.rmtree(f, ignore_errors=True)
        else:
            os.remove(f)
    for f in glob.glob(pixidir + "/envs/*"):
        shutil.rmtree(f, ignore_errors=True)

def create_pixi_dockerfile():
    with open("Dockerfile", "w") as dockerfile:
        dockerfile.write("FROM ghcr.io/prefix-dev/pixi:bookworm-slim\n\nRUN ")
        cmd = []
        for env,deps in environ.items():
            if env == "report":
                channels = ["conda-forge", "r"]
            else:
                channels = ["conda-forge", "bioconda"]
            _chan = " ".join([f"--channel {i}" for i in channels]).split()
            if env == "deconvolution":
                cmd.append(
                    " ".join(["pixi", "global", "install", *_chan, "--environment", env, "--expose", "QuickDeconvolution", *deps])
                )
            else:
                cmd.append(
                    " ".join(["pixi", "global", "install", *_chan, "--environment", env, *deps])
                )
        cmd.append("rm -rf ~/.cache/rattler")    
        dockerfile.write(' &&\\ \n\t'.join(cmd))

def create_pixi_toml():
    with open("Dockerfile", "w") as dockerfile:
        dockerfile.write("FROM ghcr.io/prefix-dev/pixi:bookworm-slim\n\n")
        dockerfile.write("COPY ./pixi.toml /root/.pixi/manifests/pixi-global.toml\n\n")
        dockerfile.write("RUN pixi global update && rm -rf ~/.cache/rattler\n\n")

    # get the name of the manifest file
    _pix = subprocess.run("pixi global list".split(), capture_output = True, text = True)
    global_manifest = _pix.stdout.splitlines()[0].split()[-1].strip("\'")
    print(global_manifest)
    # clear out the manifest
    with open(global_manifest, "w") as toml:
        toml.write("version = 1\n\n")
    reset_pixi_global()

    for env,deps in environ.items():
        if env == "report":
            channels = ["conda-forge", "r"]
        else:
            channels = ["conda-forge", "bioconda"]
        _chan = " ".join([f"--channel {i}" for i in channels]).split()
        if env == "deconvolution":
            _pix = subprocess.run(
                ["pixi", "global", "install", *_chan, "--environment", env, "--expose", "QuickDeconvolution", *deps]
            )
            if _pix.returncode > 0:
                print(_pix.stderr)
                sys.exit(1)
        else:
            _pix = subprocess.run(
                ["pixi", "global", "install", *_chan, "--environment", env, *deps]
            )
            if _pix.returncode > 0:
                print(_pix.stderr)
                sys.exit(1)

    # get the manifest file and copy to this directory
    shutil.copy(global_manifest, "resources/pixi.toml")
    # clean up global packages again
    reset_pixi_global()


