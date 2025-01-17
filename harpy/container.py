"""Command to regenerate Dockerfile for container building"""

import os
import sys
import subprocess
from pathlib import Path
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

    with open("Dockerfile.raw", "w", encoding = "utf-8") as dockerraw:
        _module = subprocess.run(
            'snakemake -s containerize.smk --containerize'.split(),
            stdout = dockerraw
        )

    with open("Dockerfile.raw", "r") as dockerraw, open("Dockerfile", "w") as dockerfile:
            # copy over the first three lines
            dockerfile.write(dockerraw.readline())
            dockerfile.write(dockerraw.readline())
            dockerfile.write(dockerraw.readline())
            dockerfile.write("\nRUN mkdir -p /conda-envs/\n")
            dockerfile.write("\nCOPY container/workflow/envs/*.yaml /conda-envs\n")
            env_hash = {}
            for line in dockerraw:
                if line.startswith("#"):
                    continue
                if line.startswith("COPY"):
                    dockercmd, env, hashname = line.split()
                    env = Path(env).stem
                    hashname = hashname.split("/")[-2]
                    env_hash[env] = hashname
            runcmds = []
            for env, _hash in env_hash.items():
                runcmds.append(f"conda env create --prefix /conda-envs/{_hash} --file /conda-envs/{env}.yaml && \\")
            runcmds.append("conda clean --all -y")
            dockerfile.write("\nRUN ")
            dockerfile.write(
                "\n\t".join(runcmds)
            )
    os.remove("containerize.smk")
