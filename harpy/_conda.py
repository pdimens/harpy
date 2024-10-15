"""Module that sets up the conda environments and dependencies for all the Harpy workflows"""

import os

def create_conda_recipes():
    """Create the YAML files of the workflow conda dependencies"""
    condachannels = ["bioconda","conda-forge"]
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
        "assembly": [
            "conda-forge::python=3"
        ],
        "metassembly": [
            "bioconda::athena_meta=1.2"
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
            "conda-forge::pandoc",
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
            "bioconda::dwgsim=1.1.14",
            "bioconda::perl-math-random",
            "bioconda::perl-inline-c",
            "bioconda::perl-parse-recdescent",
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

    # post-deployment scripts
    with open(".harpy_envs/assembly.post-deploy.sh", "w", encoding="utf-8") as shellscript:
        shellscript.write("wget -O .spades.tar.gz https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz\n")
        shellscript.write("tar -xvzf .spades.tar.gz && rm .spades.tar.gz\n")
        shellscript.write("mv SPAdes-4.0.0-Linux/bin/* ${CONDA_PREFIX}/bin && mv SPAdes-4.0.0-Linux/share/* ${CONDA_PREFIX}/share\n")
        shellscript.write("rm -r SPAdes-4.0.0-Linux\n")