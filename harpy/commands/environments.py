"""Command to regenerate Dockerfile for container building"""

import os
import shutil
import rich_click as click
from harpy.common.environments import HarpyEnvs
from harpy.common.workflow import Workflow

env_list = ['all'] + list(HarpyEnvs().environments().keys())

@click.command(hidden = True)
@click.help_option('--help', hidden = True)
@click.argument('env', required = True, type= click.Choice(env_list))
def containerize(env):
    """
    Configure the harpy container

    **INTERNAL USE ONLY**. Used to recreate all the conda environments required
    by the workflows and build a dockerfile from that.
    """
    HarpyEnvs().prepare_container(env)

@click.group()
@click.help_option('--help', hidden = True)
def deps():
    """
    Locally install workflow dependencies

    These commands are intended only for situations on HPCs where `conda` cannot
    be installed or the worker nodes do not have internet access
    to download conda/apptainer workflow dependencies.
    """

@click.command(no_args_is_help = True)
@click.help_option('--help', hidden = True)
@click.argument('workflows', required = True, type= click.Choice(env_list), nargs = -1)
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
    - impute
    - metassembly
    - phase
    - preprocess
    - qc
    - variants
    """
    workflow = Workflow("localenv", "environments.smk", "localenv/", None, False, 1)
    # if "all" was mixed with other workflows, default to just all and avoid doubling up
    _henv = HarpyEnvs()
    _henv.write_recipes(workflow.output_directory, ["all"])
    if "all" in workflows:
        workflows = list(_henv.environments().keys())
    workflow.fetch_snakefile()

    config_params = "--config"
    if "assembly" in workflows:
        config_params += " spades=True"
    config_params += " envs=[\"" + "\",\"".join(workflows) + "\"]"
    workflow.snakemake_cmd_relative = " ".join(["snakemake", "-s", os.path.join(workflow.workflow_directory, "workflow.smk"), "--sdm", "conda", "--cores 2", "--conda-prefix ../.environments", "--conda-cleanup-pkgs cache", "--directory localenv", config_params])
    workflow.launch()
    shutil.rmtree(workflow.output_directory, ignore_errors = True)

@click.command()
@click.help_option('--help', hidden = True)
def container():
    """
    Install workflow dependency containers

    Manually pull the harpy dependency container from dockerhub and convert it
    into an Apptainer .sif. To use, run this command again without arguments.
    """
    workflow = Workflow("localcontainer", "environments.smk", "localenv/", None, True, 1)
    workflow.fetch_snakefile()
    workflow.snakemake_cmd_relative = " ".join(["snakemake", "-s", os.path.join(workflow.workflow_directory, "workflow.smk"), "--sdm", "apptainer", "--cores 2", "--apptainer-prefix ../.environments", "--directory localenv"])
    workflow.launch()
    shutil.rmtree(workflow.output_directory, ignore_errors = True)

deps.add_command(conda)
deps.add_command(container)