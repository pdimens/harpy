"""Creates environment recipes for all the Harpy workflows"""

import os
import sys
import yaml
from rich import box
from rich.table import Table
from .printing import print_error

def create_conda_recipes(outdir: str, envs: list= []) -> None:
    """Create the YAML files of the workflow conda dependencies"""
    environ = {
        "align" : [
            "bioconda::bwa",
            "bioconda::ema",
            "bioconda::samtools=1.22",
            "bioconda::seqtk",
            "bioconda::strobealign",
            "bioconda::tabix",
            "conda-forge::icu",
            "conda-forge::libzlib",
            "conda-forge::xz"
        ],
        "assembly" : [
            "bioconda::arcs",
            "bioconda::bwa",
            "bioconda::cloudspades",
            "bioconda::links",
            "bioconda::quast",
            "bioconda::busco",
            "bioconda::samtools",
            "bioconda::tigmint"
        ],
        "deconvolution" : [
            "bioconda::quickdeconvolution"
        ],

        "demultiplex": [
            "bioconda::dmox>=0.2"
        ],
        "metassembly": [
            "bioconda::athena_meta=1.2"
        ],
        "phase" : [
            "bioconda::hapcut2",
            "bioconda::whatshap"
        ],
        "qc" : [
            "conda-forge::click=8.2.1",
            "bioconda::falco=1.2.5",
            "bioconda::fastp",
            "bioconda::multiqc=1.30",
            "bioconda::pysam=0.23"
        ],
        "report" : [
            "conda-forge::quarto",
            "conda-forge::r-dt",
            "conda-forge::r-dplyr",
            "conda-forge::r-highcharter",
            "conda-forge::r-magrittr",
            "conda-forge::r-plotly",
            "conda-forge::r-scales",
            "conda-forge::r-tidyr",
            "conda-forge::r-viridislite", 
            "conda-forge::r-xml2",
            "r::r-biocircos"
        ],
        "simulations" : [
            "bioconda::mimick>=2.3",
            "bioconda::simug>1.0.0",
        ],
        "spades" : [
            "conda-forge::python=3"
        ],
        "stitch" : [
            "bioconda::r-stitch>=1.8"
        ],
        "variants" : [
            "bioconda::bcftools=1.22",
            "bioconda::freebayes=1.3.9",
            "bioconda::leviathan",
            "bioconda::naibr-plus",
            "conda-forge::setuptools"
        ]
    }
    _out = os.path.join(outdir, "workflow", "envs")
    os.makedirs(_out, exist_ok = True)
    # if none provided, use all
    if not envs:
        envs = environ.keys()

    for i in envs:
        try:
            env_dict = {
                "name" : i,
                "channels" : ["conda-forge", "bioconda"],
                "dependencies": environ[i]
            }
            if i == "r":
                env_dict["channels"].append("r")
        except KeyError:
            sys.stderr.write(f"Key '{i}' is not an available conda environment name. The options are: " + ", ".join(environ.keys()))
            sys.exit(1)
        with open(os.path.join(_out, f"{i}.yaml"), "w", encoding="utf-8") as recipe:
            yaml.dump(env_dict, recipe, default_flow_style= False, sort_keys=False, width=float('inf'), indent=2)

    if "spades" in envs:
        # post-deployment script
        with open(os.path.join(_out, "spades.post-deploy.sh"), "w", encoding="utf-8") as shellscript:
            shellscript.write("wget -O .spades.tar.gz https://github.com/ablab/spades/releases/download/v4.1.0/SPAdes-4.1.0-Linux.tar.gz\n")
            shellscript.write("tar -xvzf .spades.tar.gz && rm .spades.tar.gz\n")
            shellscript.write("mv SPAdes-4.1.0-Linux/bin/* ${CONDA_PREFIX}/bin && mv SPAdes-4.1.0-Linux/share/* ${CONDA_PREFIX}/share\n")
            shellscript.write("rm -r SPAdes-4.1.0-Linux\n")

def check_environments(dirpath: str, envs: list) -> None:
    """Check that the provided dir exists and contains the necessary environment definitions"""
    if not os.path.exists(f"{dirpath}/workflow/envs"):
        print_error("missing conda files", "This working directory does not contain the expected directory of conda environment definitions ([blue bold]workflow/envs/[/])\n  - use [green bold]--conda[/] to recreate it")
    envlist = os.listdir(f"{dirpath}/workflow/envs")
    errcount = 0
    errtable = Table(show_footer=True, box=box.SIMPLE)
    errtable.add_column("File", justify="left", no_wrap=True)
    errtable.add_column("Status", justify="center")
    for i in envs:
        if f"{i}.yaml" in envlist:
            errtable.add_row(f"[dim]{i}.yaml", "[dim]present")
        else:
            errcount += 1
            errtable.add_row(f"[yellow bold]{i}.yaml", "[yellow bold]missing")
    if errcount > 0:
        print_error(
            "missing environment files",
            f"The directory [blue]{dirpath}/workflows/envs[/] is missing [yellow bold]{errcount}[/] of the expected conda environment definition files.",
            "Check that the names conform to Harpy's expectations, otherwise you can recreate this directory using the [green bold]--conda[/] option.",
            "Expected environment files",
            errtable
            )
