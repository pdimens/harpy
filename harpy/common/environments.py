"""Creates environment recipes for all the Harpy workflows"""

import copy
import os
import shutil
import subprocess
import sys
import yaml
from rich import box
from rich.table import Table
from .printing import HarpyPrint
from harpy.common.version import VERSION

class HarpyEnvs():
    '''The class that holds conda and pixi environments and means to create container dockerfiles'''
    def __init__(self):
        self.__environments__: dict = {
        "align" : [
            "bioconda::bwa-mem2",
            "bioconda::bwa",
            "bioconda::samtools=1.22",
            "bioconda::seqtk",
            "bioconda::strobealign",
            "bioconda::tabix"
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
        "metassembly": [
            "bioconda::athena_meta=1.2"
        ],
        "phase" : [
            "bioconda::hapcut2",
            "bioconda::whatshap"
        ],
        "qc" : [
            "conda-forge::click=8.2.1",
            "bioconda::cutadapt",
            "bioconda::dmox>=0.2",
            "bioconda::falco=1.2.5",
            "bioconda::fastp",
            "bioconda::multiqc=1.30",
            "bioconda::pysam=0.23",
            "bioconda::quickdeconvolution"
        ],
        "impute" : [
            "bioconda::r-stitch>=1.8.4"
        ],
        "variants" : [
            "bioconda::bcftools=1.22",
            "bioconda::freebayes=1.3.9",
            "bioconda::leviathan",
            "bioconda::naibr-plus",
            "conda-forge::setuptools"
        ]
    }

        self.dockerfile: str = """
FROM ghcr.io/prefix-dev/pixi:0.62.0 AS build

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

    def environments(self, is_pixi: bool = False) -> dict:
        '''
        Return the `dict` of software environments. Use `is_pixi=True` to remove the
        `channel::` prefix from software.
        '''
        d = copy.copy(self.__environments__)
        if is_pixi:
            for i in d:
                d[i] = [j.split("::")[-1] for j in d[i]]
        return d

    def write_recipes(self, outdir: str, envs: list= []) -> None:
        """Create the YAML files of the workflow conda dependencies"""
        environ = self.__environments__
        _out = os.path.join(outdir, "workflow", "envs")
        os.makedirs(_out, exist_ok = True)
        # if none provided, use all
        if not envs:
            envs = list(environ.keys())

        for i in envs:
            try:
                env_dict = {
                    "name" : i,
                    "channels" : ["conda-forge", "bioconda"],
                    "dependencies": environ[i]
                }
                if i == "report":
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

    def prepare_container(self, env):
        '''
        Using the defined environments, create a folder (or series of folders) with a dockerfile
        and pixi.toml file to create one of the environments.
        '''
        environ = self.environments(True)

        if env != "all":
            environ = {env: environ.get(env)}

        for env,deps in environ.items():
            shutil.rmtree(f"container/{env}", ignore_errors=True)
            os.makedirs(f"container/{env}", exist_ok=True)
            with open(f"container/{env}/Dockerfile", "w") as dockerfile:
                dockerfile.write(self.dockerfile)
            if env == "report":
                subprocess.run(
                    f"pixi init container/{env} -c conda-forge -c r".split(),
                    check = True
                )
            else:
                subprocess.run(
                    f"pixi init container/{env} -c conda-forge -c bioconda".split(),
                    check = True    
                )

            with open(f"container/{env}/pixi.toml", "r") as toml:
                with open(f"container/{env}/pixi.fix.toml", "w") as out:
                    for line in toml:
                        if line.startswith("version"):
                            line = f"version = \"{VERSION}\"\n"
                        out.write(line)

            os.remove(f"container/{env}/pixi.toml")
            shutil.move(f"container/{env}/pixi.fix.toml", f"container/{env}/pixi.toml")

            subprocess.run(
                ["pixi", "add", "--no-progress", "--manifest-path", f"container/{env}/pixi.toml", *deps],
                check = True
            )

            shutil.rmtree("container/.pixi", ignore_errors=True)


def check_environments(dirpath: str, envs: list) -> None:
    """Check that the provided dir exists and contains the necessary environment definitions"""
    if not os.path.exists(f"{dirpath}/workflow/envs"):
        HarpyPrint().error("missing conda files", "This working directory does not contain the expected directory of conda environment definitions ([blue bold]workflow/envs/[/])\n  - use [green bold]--conda[/] to recreate it")
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
        HarpyPrint().error(
            "missing environment files",
            f"The directory [blue]{dirpath}/workflows/envs[/] is missing [yellow bold]{errcount}[/] of the expected conda environment definition files.",
            "Check that the names conform to Harpy's expectations, otherwise you can recreate this directory using the [green bold]--conda[/] option.",
            "Expected environment files",
            errtable
            )
