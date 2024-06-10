"""Module to bypass Harpy and run snakemake"""

import os
import sys
import subprocess
import yaml
import rich_click as click
from rich.markdown import Markdown
from .printfunctions import print_error, print_onstart
from .conda_deps import generate_conda_deps

#docstring = {
#    "harpy run": [
#        {
#            "name": "Submission Systems",
#            "commands": ["slurm", "htcondor", "lsf", "googlebatch", "generic"],
#        },
#    ]
#}

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/other")
@click.option('-c', '--conda',  is_flag = True, default = False, help = 'Recreate the conda environments into .harpy_envs/')
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False), nargs=1)
def run(directory, conda):
    """
    Run a workflow from an existing Harpy directory

    In the event you need to run the Snakemake workflow present in a Harpy output directory
    (e.g. `Align/bwa`) without Harpy rewriting any of the configuration files, this command
    bypasses all the initialization steps of Harpy workflows and executes the Snakemake command
    present in `directory/workflow/config.yaml`. The only requirements are:
    1. the target directory has `workflow/config.yaml` present in it
    2. the current working directory has a `.harpy_envs` directory with the necessary environments defined
    """
    errtext = []
    directory = directory.rstrip("/")
    if not os.path.exists(f"{directory}/workflow/config.yaml"):
        errtext.append(f"- Target directory `{directory}` does not contain `workflow/config.yaml`")
    if not os.path.exists(".harpy_envs"):
        errtext.append("- Current directory does not contain `.harpy_envs/`\n  - use `--conda` to recreate it")
    else:
        if len(os.listdir(".harpy_envs")) == 0:
            errtext.append("- The directory of conda environments (`.harpy_envs/`) is empty`\n  - use `--conda` to recreate it")
    if errtext:
        print_error(Markdown("\n".join(errtext)))
        sys.exit(1)
    
    with open(f"{directory}/workflow/config.yaml", 'r') as f:
        harpy_config = yaml.full_load(f)
        command = harpy_config["workflow_call"]
    
    if conda:
        generate_conda_deps()
    
    print_onstart(f"Output Directory: {directory}", harpy_config["workflow"])
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)