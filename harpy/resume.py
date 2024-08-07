"""Module to bypass Harpy and run snakemake"""

import os
import sys
import subprocess
import yaml
import rich_click as click
from .validations import check_envdir
from .printfunctions import print_error, print_onstart
from .conda_deps import generate_conda_deps

@click.command(no_args_is_help = True, epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/other")
@click.option('-c', '--conda',  is_flag = True, default = False, help = 'Recreate the conda environments into .harpy_envs/')
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False, readable=True), nargs=1)
def resume(directory, conda):
    """
    Resume a workflow from an existing Harpy directory

    In the event you need to run the Snakemake workflow present in a Harpy output directory
    (e.g. `Align/bwa`) without Harpy rewriting any of the configuration files, this command
    bypasses all the initialization steps of Harpy workflows and executes the Snakemake command
    present in `directory/workflow/config.yaml`. It will reuse an existing `.harpy_envs/` folder,
    otherwise use `--conda` to create one.

    The only requirements is:
    - the target directory has `workflow/config.yaml` present in it
    """
    directory = directory.rstrip("/")
    if not os.path.exists(f"{directory}/workflow/config.yaml"):
        print_error(f"Target directory [blue bold]{directory}[/blue bold] does not contain [blue bold]workflow/config.yaml[/blue bold]")
        sys.exit(1)
    if conda:
        generate_conda_deps()
    else:
        check_envdir(".harpy_envs")
    with open(f"{directory}/workflow/config.yaml", 'r', encoding="utf-8") as f:
        harpy_config = yaml.full_load(f)
        command = harpy_config["workflow_call"]
    
    print_onstart(f"Output Directory: {directory}", "resume: " + harpy_config["workflow"])
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)



