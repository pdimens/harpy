"""Module to bypass Harpy and run snakemake"""

import os
import re
import sys
import yaml
import rich_click as click
from ._conda import check_environments
from ._cli_types_generic import convert_to_int
from ._printing import print_error, workflow_info
from ._launch import launch_snakemake
from ._misc import snakemake_log, write_workflow_config
from ._conda import create_conda_recipes

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/other")
@click.option('-c', '--conda',  is_flag = True, default = False, help = 'Recreate the conda environments')
@click.option('-r', '--relative',  is_flag = True, default = False, help = 'Call Snakemake with relative paths')
@click.option('-t', '--threads', type = click.IntRange(2, 999, clamp = True), help = 'Change the number of threads (>1)')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False, readable=True, resolve_path=True), nargs=1)
def resume(directory, conda, relative, threads, quiet):
    """
    Resume a Harpy workflow from an existing directory

    In the event you need to run the Snakemake workflow present in a Harpy output directory
    (e.g. `Align/bwa`) without Harpy redoing validations and rewriting any of the configuration files,
    this command bypasses all the preprocessing steps of Harpy workflows and executes the Snakemake command
    present in `directory/workflow/config.yaml`. It will reuse an existing `workflow/envs/` folder
    for conda environments, otherwise use `--conda` to create one.

    The only requirements are:
    - the target directory has `workflow/config.yaml` present in it
    - the targest directory has `workflow/envs/*.yaml` present in it
    """
    if not os.path.exists(f"{directory}/workflow/config.yaml"):
        print_error("missing snakemake config", f"Target directory [blue]{directory}[/] does not contain the file [bold]workflow/config.yaml[/]")
        sys.exit(1)
    if not os.path.exists(f"{directory}/workflow/config.harpy.yaml"):
        print_error("missing workflow config", f"Target directory [blue]{directory}[/] does not contain the file [bold]workflow/config.harpy.yaml[/]")
        sys.exit(1)
    
    with open(f"{directory}/workflow/config.harpy.yaml", 'r', encoding="utf-8") as f:
        harpy_config = yaml.full_load(f)
    conda_envs = harpy_config["conda_environments"] 
    if conda:
        create_conda_recipes(directory, conda_envs)
    else:
        check_environments(directory, conda_envs)
    
    workflow = harpy_config["workflow"]
    sm_log = snakemake_log(directory, workflow)
    if sm_log != harpy_config["snakemake"]["log"]:
        harpy_config["snakemake"]["log"] = sm_log
        ## overwrite config file with this updated one, where only the snakemake log has been changed ##
        write_workflow_config(harpy_config, directory)

    command = harpy_config["snakemake"]["relative"]
    if threads:
        command = re.sub(r"--cores \d+", f"--cores {threads}", command)
    start_text = workflow_info(
        ("Workflow:", workflow.replace("_", " ")),
        ("Output Folder:", directory + "/"),
        ("Workflow Log:", sm_log.replace(f"{directory}/", "") + "[dim].gz")
    )
    launch_snakemake(command_rel, workflow, start_text, directory, sm_log, quiet, f'workflow/{workflow}.summary')



