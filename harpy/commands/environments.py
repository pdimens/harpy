"""Command to regenerate Dockerfile for container building"""

import os
import shutil
import subprocess
import rich_click as click
from harpy.common.conda import create_conda_recipes
from harpy.common.workflow import Workflow

@click.command(hidden = True)
def containerize():
    """
    Configure conda and docker environments

    **INTERNAL USE ONLY**. Used to recreate all the conda environments required
    by the workflows and build a dockerfile from that.
    """
    workflow = Workflow("container", "environments.smk", "container", 1)
    workflow.fetch_snakefile()
    create_conda_recipes("container")
    
    with open("Dockerfile", "w", encoding = "utf-8") as dockerraw:
        _module = subprocess.run(
            'snakemake -s container/workflow/workflow.smk --containerize --directory container'.split(),
            stdout = dockerraw
        )

    #with open("Dockerfile.raw", "r") as dockerraw, open("Dockerfile", "w") as dockerfile:
    #        # copy over the first three lines
    #        dockerfile.write(dockerraw.readline())
    #        dockerfile.write(dockerraw.readline())
    #        dockerfile.write(dockerraw.readline())
    #        dockerfile.write("\nCOPY container/workflow/envs/*.yaml /\n")
    #        env_hash = {}
    #        for line in dockerraw:
    #            if line.startswith("#"):
    #                continue
    #            if line.startswith("COPY"):
    #                dockercmd, env, hashname = line.split()
    #                env = Path(env).stem
    #                hashname = hashname.split("/")[-2]
    #                env_hash[env] = hashname
    #        runcmds = []
    #        for env, _hash in env_hash.items():
    #            runcmds.append(f"conda env create --prefix /conda-envs/{_hash} --file /{env}.yaml && \\")
    #        runcmds.append("conda clean --all -y")
    #        dockerfile.write("\nRUN ")
    #        dockerfile.write(
    #            "\n\t".join(runcmds)
    #        )
    #os.remove("Dockerfile.raw")

@click.group(options_metavar='')
def deps():
    """
    Locally install workflow dependencies

    These commands are intended only for situations on HPCs where `conda` cannot
    be installed or the worker nodes do not have internet access
    to download conda/apptainer workflow dependencies.
    """

@click.command(no_args_is_help = True)
@click.argument('workflows', required = True, type= click.Choice(["all", "align", "assembly", "metassembly", "phase", "qc", "report", "simulations", "stitch", "variants"]), nargs = -1)
def conda(workflows):
    """
    Install workflow dependencies via conda

    Create the conda environments required by Harpy's workflows (e.g. `phase`).
    ONLY use this for specific HPC configurations where worker nodes
    do not have internet access to let snakemake install conda packages. 
    Provide any combination of: 
    - all
    - align
    - assembly
    - metassembly
    - phase
    - qc
    - r
    - simulations
    - stitch
    - variants
    """
    workflow = Workflow("localenv", "environments.smk", "localenv/", 1)
    # if "all" was mixed with other workflows, default to just all and avoid doubling up
    create_conda_recipes(workflow.output_directory)
    if "all" in workflows:
        workflows = ["align", "assembly", "metassembly", "phase", "qc", "r", "simulations", "stitch", "variants"] 
    workflow.fetch_snakefile()

    config_params = "--config"
    if "assembly" in workflows:
        config_params += " spades=True"
    config_params += " envs=[\"" + "\",\"".join(workflows) + "\"]"
    workflow.snakemake_cmd_relative = " ".join(["snakemake", "-s", os.path.join(workflow.workflow_directory, "workflow.smk"), "--sdm", "conda", "--cores 2", "--conda-prefix ../.environments", "--conda-cleanup-pkgs cache", "--directory localenv", config_params])
    workflow.launch()
    shutil.rmtree(workflow.output_directory, ignore_errors = True)

@click.command(context_settings={"help_option_names" : ["-h", "--help"]})
def container():
    """
    Install workflow dependency container

    Manually pull the harpy dependency container from dockerhub and convert it
    into an Apptainer .sif. To use, run this command again without arguments.
    """
    workflow = Workflow("localcontainer", "environments.smk", "localenv/", 1)
    workflow.fetch_snakefile()
    workflow.snakemake_cmd_relative = " ".join(["snakemake", "-s", os.path.join(workflow.workflow_directory, "workflow.smk"), "--sdm", "conda apptainer", "--cores 2", "--apptainer-prefix ../.environments", "--directory localenv"])
    workflow.launch()
    shutil.rmtree(workflow.output_directory, ignore_errors = True)

deps.add_command(conda)
deps.add_command(container)