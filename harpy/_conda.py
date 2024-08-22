"""Module that sets up the conda environments and dependencies for all the Harpy workflows"""

import os
import sys
import yaml
from rich import box
from rich.table import Table
from rich import print as rprint
from ._printing import print_error, print_solution_with_culprits

def create_conda_recipes(outdir: str, envs: list=None) -> None:
    """Create the YAML files of the workflow conda dependencies"""
    condachannels = ["bioconda","conda-forge","defaults"]
    environ = {
        "align" : [
            "bioconda::bwa",
            "bioconda::ema",
            "bioconda::samtools=1.20",
            "bioconda::seqtk",
            "bioconda::strobealign",
            "bioconda::tabix",
            "conda-forge::icu",
            "conda-forge::libzlib",
            "conda-forge::xz"
            ],
        "metassembly" : [
            "bioconda::spades=4.0"
        ],
        "athena": [
            "bioconda::athena_meta"
        ],
        "phase" : [
            "bioconda::hapcut2",
            "bioconda::whatshap"
            ],
        "qc" : [
            "bioconda::falco",
            "bioconda::fastp",
            "bioconda::multiqc=1.22",
            "bioconda::pysam=0.22",
            "bioconda::quickdeconvolution"
            ],
        "r" : [
            "conda-forge::r-dt",
            "conda-forge::r-dplyr",
            "conda-forge::r-flexdashboard",
            "conda-forge::r-ggplot2",
            "conda-forge::r-highcharter",
            "conda-forge::r-magrittr",
            "conda-forge::r-plotly",
            "conda-forge::r-scales",
            "conda-forge::r-stringi",
            "conda-forge::r-tidyr",
            "conda-forge::r-viridislite", 
            "conda-forge::r-xml2",
            "r::r-biocircos"
            ],
        "simulations" : [
            "alienzj::msort",
            "bioconda::dwgsim",
            "bioconda::perl-math-random",
            "bioconda::perl-inline-c",
            "bioconda::perl-parse-recdescent",
            "bioconda::simug>1.0.0",
            "conda-forge::numpy",
            "conda-forge::perl"
            ],
        "variants": [
            "bioconda::bcftools=1.20",
            "bioconda::freebayes=1.3.6",
            "bioconda::leviathan",
            "bioconda::naibr-plus"
            ],
        "stitch" : [
            "bioconda::r-stitch=1.6.10"
            ]
    }

    os.makedirs(".harpy_envs", exist_ok = True)
    # overwrites existing
    for env,deps in environ.items():
        with open(f".harpy_envs/{env}.yaml", mode="w", encoding="utf-8") as yml:
            yml.write(f"name: {env}\n")
            yml.write("channels:\n  - ")
            yml.write("\n  - ".join(condachannels))
            yml.write("\ndependencies:\n  - ")
            yml.write("\n  - ".join(deps) + "\n")