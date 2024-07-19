"""Module that sets up the conda environments and dependencies for all the Harpy workflows"""

import os

def generate_conda_deps():
    """Create the YAML files of the workflow conda dependencies"""
    condachannels = ["bioconda","conda-forge","defaults"]
    environ = {
        "align": [
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
            "conda-forge::r-flexdashboard",
            "conda-forge::r-ggplot2",
            "conda-forge::r-highcharter",
            "conda-forge::r-plotly",
            "conda-forge::r-tidyr",
            "conda-forge::r-xml2",
            "r::r-biocircos"
            ],
        "simulations" : [
            "alienzj::msort",
            "bioconda::dwgsim",
            "bioconda::perl-math-random",
            "bioconda::perl-inline-c",
            "bioconda::perl-parse-recdescent",
            "conda-forge::numpy",
            "conda-forge::perl"
            ],
        "snp": [
            "bioconda::bcftools=1.20",
            "bioconda::freebayes=1.3.6"
            ],
        "stitch" : [
            "bioconda::r-stitch=1.6.10"
            ],
        "sv": [
            "bioconda::leviathan",
            "bioconda::naibr-plus"
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
