"""Command to regenerate Dockerfile for container building"""

import os
import sys
import shutil
import subprocess
from pathlib import Path
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, snakemake_log

@click.command(hidden = True)
def containerize():
    """
    Configure conda and docker environments

    **INTERNAL USE ONLY**. Used to recreate all the conda environments required
    by the workflows and build a dockerfile from that.
    """
    create_conda_recipes("container")
    fetch_rule("container/workflow", "environments.smk")

    with open("Dockerfile.raw", "w", encoding = "utf-8") as dockerraw:
        _module = subprocess.run(
            'snakemake -s workflow/workflow.smk --containerize --directory container'.split(),
            stdout = dockerraw
        )

    with open("Dockerfile.raw", "r") as dockerraw, open("Dockerfile", "w") as dockerfile:
            # copy over the first three lines
            dockerfile.write(dockerraw.readline())
            dockerfile.write(dockerraw.readline())
            dockerfile.write(dockerraw.readline())
            #dockerfile.write("\nRUN mkdir -p /conda-envs/\n")
            dockerfile.write("\nCOPY container/workflow/envs/*.yaml /\n")
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
                runcmds.append(f"conda env create --prefix /conda-envs/{_hash} --file /{env}.yaml && \\")
            runcmds.append("conda clean --all -y")
            dockerfile.write("\nRUN ")
            dockerfile.write(
                "\n\t".join(runcmds)
            )
    os.remove("Dockerfile.raw")
    os.remove("environments.smk")

@click.command(hidden = True)
@click.argument('workflows', required = True, type= click.Choice(["all", "align", "assembly", "metassembly", "phase", "qc", "r", "simulations", "stitch", "variants"]), nargs = -1)
def localenv(workflows):
    """
    Install Harpy workflow dependencies via conda

    **INTERNAL USE ONLY**. Used to recreate the conda environments required
    by the workflows. Provide any combination of: 
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
    output_dir = "localenv/"
    sm_log = snakemake_log(output_dir, "localenv")
    # if "all" was mixed with other workflows, default to just all and avoid doubling up
    if "all" in workflows:
        create_conda_recipes(output_dir, workflows)
    else:
        create_conda_recipes(output_dir)
    fetch_rule(os.path.join(output_dir, 'workflow'), "environments.smk")
    command = " ".join(["snakemake", "-s", os.path.join(output_dir, "workflow", "workflow.smk"), "--sdm", "conda", "--cores 2", "--conda-prefix .environments", "--conda-cleanup-pkgs cache", "--directory .", "--config spades=True"])
    launch_snakemake(command, "localenv", "", output_dir, sm_log, 1, "workflow/localenv.summary")
    shutil.rmtree(output_dir, ignore_errors = True)
