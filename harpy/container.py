"""Command to regenerate Dockerfile for container building"""

import os
import sys
import subprocess
import rich_click as click
from ._conda import create_conda_recipes
from ._misc import fetch_rule

@click.command(no_args_is_help = False, hidden = True)
def containerize():
    """
    Configure conda and docker environments

    **INTERNAL USE ONLY**. Used to recreate all the conda environments required
    by the workflows and build a dockerfile from that.
    """
    create_conda_recipes("container")
    fetch_rule(os.getcwd(), "containerize.smk")

    with open("Dockerfile", "w", encoding = "utf-8") as dockerfile:
        _module = subprocess.run(
            'snakemake -s containerize.smk --containerize --conda-cleanup-pkgs cache'.split(),
            stdout = dockerfile
        )
    os.remove("containerize.smk")
    sys.exit(_module.returncode)
