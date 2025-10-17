#! /usr/bin/env python

import os
import shutil
import subprocess

def create_pixi_toml():
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
            "spades",
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
            "mimick>=2.3",
            "simug>=1.0.1",
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

    for env,deps in environ.items():
        if env == "report":
            channels = ["conda-forge", "r"]
        else:
            channels = ["conda-forge", "bioconda"]
        _chan = " ".join([f"--channel {i}" for i in channels]).split()
        
        if env == "deconvolve":
            subprocess.run(["pixi", "global", "install", *_chan, "--environment", env, "--expose", "quickdeconvolve", *deps])
        else:
            subprocess.run(["pixi", "global", "install", *_chan, "--environment", env, *deps])
    shutil.copy2(os.path.expanduser("~/.pixi/manifests/pixi-global.toml"), "resources/pixi.toml")
