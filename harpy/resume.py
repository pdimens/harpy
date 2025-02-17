"""Module to bypass Harpy and run snakemake"""

import os
import re
import sys
import yaml
import rich_click as click
from ._conda import check_environments
from ._printing import print_error, workflow_info
from ._launch import launch_snakemake
from ._misc import snakemake_log
from ._conda import create_conda_recipes

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/other")
@click.option('-c', '--conda',  is_flag = True, default = False, help = 'Recreate the conda environments')
@click.option('-t', '--threads', type = click.IntRange(min = 2, max_open = True), help = 'Change the number of threads (>1)')
@click.option('--quiet',  is_flag = True, default = False, help = 'Don\'t show output text while running')
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False, readable=True), nargs=1)
def resume(directory, conda, threads, quiet):
    """
    Resume a workflow from an existing Harpy directory

    In the event you need to run the Snakemake workflow present in a Harpy output directory
    (e.g. `Align/bwa`) without Harpy rewriting any of the configuration files, this command
    bypasses all the preprocessing steps of Harpy workflows and executes the Snakemake command
    present in `directory/workflow/config.yaml`. It will reuse an existing `workflow/envs/` folder
    for conda environments, otherwise use `--conda` to create one.

    The only requirements are:
    - the target directory has `workflow/config.yaml` present in it
    - the targest directory has `workflow/envs/*.yaml` present in it
    """
    directory = directory.rstrip("/")
    if not os.path.exists(f"{directory}/workflow/config.yaml"):
        print_error("missing config file", f"Target directory [blue]{directory}[/blue] does not contain the file [bold]workflow/config.yaml[/bold]")
        sys.exit(1)
    with open(f"{directory}/workflow/config.yaml", 'r', encoding="utf-8") as f:
        harpy_config = yaml.full_load(f)
    conda_envs = harpy_config["conda_environments"] 
    if conda:
        create_conda_recipes(directory, conda_envs)
    else:
        check_environments(directory, conda_envs)
    
    workflow = harpy_config["workflow"].replace(" ", "_")
    sm_log = snakemake_log(directory, workflow)
    command = harpy_config["workflow_call"] + f" --config snakemake_log={sm_log}"
    if threads:
        command = re.sub(r"--threads \d+", f"--threads {threads}", command)
    start_text = workflow_info(
        ("Output Folder:", directory + "/"),
        ("Workflow Log:", sm_log.replace(f"{directory}/", "") + "[dim].gz")
    )
    launch_snakemake(command, workflow, start_text, directory, sm_log, int(quiet), f'workflow/{workflow.replace("_", ".")}.summary')



