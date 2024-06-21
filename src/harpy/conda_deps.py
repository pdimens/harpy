"""Module that sets up the conda environments and dependencies for all the Harpy workflows"""

import os

def generate_conda_deps():
    """Create the YAML files of the workflow conda dependencies"""
    condachannels = ["bioconda","conda-forge","defaults"]
    environ = {
        "qc" : ["bioconda::falco", "bioconda::fastp", "bioconda::multiqc", "bioconda::pysam=0.22"],
        "align": ["bioconda::bwa", "bioconda::ema","bioconda::strobealign", "conda-forge::icu","conda-forge::libzlib", "bioconda::samtools=1.20", "bioconda::seqtk", "bioconda::tabix", "conda-forge::xz"],
        "snp": ["bioconda::bcftools=1.20", "bioconda::freebayes=1.3.6"],
        "sv": ["bioconda::leviathan", "bioconda::naibr-plus"],
        "phase" : ["bioconda::hapcut2", "bioconda::whatshap"],
        "simulations" : ["conda-forge::perl", "bioconda::perl-math-random", "bioconda::perl-inline-c", "bioconda::perl-parse-recdescent", "conda-forge::numpy", "bioconda::dwgsim", "alienzj::msort"],
        "r" : ["conda-forge::r-xml2", "conda-forge::r-highcharter", "conda-forge::r-circlize", "r::r-biocircos", "conda-forge::r-dt", "conda-forge::r-flexdashboard", "conda-forge::r-ggplot2", "conda-forge::r-ggridges", "conda-forge::r-plotly", "conda-forge::r-tidyr", "bioconda::r-stitch"]
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
