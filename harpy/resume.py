"""Module to bypass Harpy and run snakemake"""

import os
import sys
import yaml
from rich import box
from rich.table import Table
import rich_click as click
from ._validations import check_envdir
from ._printing import print_error
from ._launch import launch_snakemake
from ._misc import snakemake_log
from ._conda import create_conda_recipes

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/other")
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
        print_error("config.yaml missing", f"Target directory [blue bold]{directory}[/blue bold] does not contain [blue bold]workflow/config.yaml[/blue bold]")
        sys.exit(1)
    if conda:
        create_conda_recipes()
    else:
        check_envdir(".harpy_envs")
    with open(f"{directory}/workflow/config.yaml", 'r', encoding="utf-8") as f:
        harpy_config = yaml.full_load(f)
    
    workflow = harpy_config["workflow"].replace(" ", "_")
    sm_log = snakemake_log(directory, workflow)
    command = harpy_config["workflow_call"] + f" --config snakemake_log={sm_log}"
    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column(header="value", justify="left")
    start_text.add_row("Output Folder:", directory + "/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{directory}/", ""))
    launch_snakemake(command, workflow, start_text, directory, sm_log, quiet)



