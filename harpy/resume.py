"""Module to bypass Harpy and run snakemake"""

import os
import sys
import yaml
import rich_click as click
from .validations import check_envdir
from .printfunctions import print_error
from .helperfunctions import snakemake_log, launch_snakemake
from .conda_deps import generate_conda_deps

@click.command(no_args_is_help = True, epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/other")
@click.option('-c', '--conda',  is_flag = True, default = False, help = 'Recreate the conda environments into .harpy_envs/')
@click.option('--quiet',  is_flag = True, default = False, help = 'Don\'t show output text while running')
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False, readable=True), nargs=1)
def resume(directory, conda, quiet):
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
    
    workflow = harpy_config["workflow"].replace(" ", "_")
    sm_log = snakemake_log(directory, workflow)
    command = harpy_config["workflow_call"] + f" --config snakemake_log={sm_log}"
    start_text = f"Output Directory: {directory}\nLog: {sm_log}"
    launch_snakemake(command, workflow, start_text, directory, sm_log, quiet)


