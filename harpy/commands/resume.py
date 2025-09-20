"""Module to bypass Harpy and run snakemake"""

from datetime import datetime
import os
import re
import yaml
import rich_click as click
from harpy.common.conda import check_environments, create_conda_recipes
from harpy.common.printing import print_error, workflow_info
from harpy.common.workflow import Workflow

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/other")
@click.option('-c', '--conda',  is_flag = True, default = False, help = 'Recreate the conda environments')
@click.option('-a', '--absolute',  is_flag = True, default = False, help = 'Call Snakemake with absolute paths')
@click.option('-t', '--threads', type = click.IntRange(2, 999, clamp = True), help = 'Change the number of threads (>1)')
@click.option('--quiet', default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False, readable=True, resolve_path=True), nargs=1)
def resume(directory, conda, absolute, threads, quiet):
    """
    Continue an incomplete Harpy workflow

    In the event you need to run the Snakemake workflow present in a Harpy output directory
    (e.g. `Align/bwa`) without Harpy redoing validations and rewriting any of the configuration files,
    this command bypasses all the preprocessing steps of Harpy workflows and executes the Snakemake command
    present in `directory/workflow/workflow.yaml`. It will reuse an existing `workflow/envs/` folder
    to validate software dependencies, otherwise use `--conda` to create a populated one.

    The only requirements are:
    - the target directory has `workflow/config.yaml` present in it
    - the targest directory has `workflow/envs/*.yaml` present in it
    """
    CONFIG_FILE = os.path.join(directory, "workflow", "workflow.yaml")
    PROFILE_FILE = os.path.join(directory, "workflow", "config.yaml")
    if not os.path.exists(PROFILE_FILE):
        print_error("missing snakemake config", f"Target directory [yellow]{directory}[/] does not contain the file [blue]workflow/config.yaml[/]")
    if not os.path.exists(CONFIG_FILE):
        print_error("missing workflow config", f"Target directory [yellow]{directory}[/] does not contain the file [blue]workflow/workflow.yaml[/]")
    
    with open(CONFIG_FILE, 'r', encoding="utf-8") as f:
        harpy_config = yaml.full_load(f)
    with open(PROFILE_FILE, 'r', encoding="utf-8") as f:
        snakemake_config = yaml.full_load(f)

    workflow = Workflow(harpy_config["workflow"], "NA", snakemake_config["directory"], quiet)
    workflow.conda = harpy_config["conda_environments"] 

    if conda:
        create_conda_recipes(directory, workflow.conda)
    else:
        check_environments(directory, workflow.conda)
    
    sm_log = os.path.join(directory, harpy_config["snakemake"]["log"])
    if os.path.exists(sm_log) or os.path.exists(sm_log + ".gz"):
        timestamp = datetime.now().strftime("%d_%m_%Y") + ".log"
        split_log = sm_log.split(".")
        _basename = ".".join(split_log[0:-3])
        incremenent = int(split_log[-3]) + 1
        harpy_config["snakemake"]["log"] = f"{_basename}.{incremenent}.{timestamp}"

    if threads:
        harpy_config["snakemake"]["absolute"] = re.sub(r"--cores \d+", f"--cores {threads}", harpy_config["snakemake"]["absolute"])
        harpy_config["snakemake"]["relative"] = re.sub(r"--cores \d+", f"--cores {threads}", harpy_config["snakemake"]["relative"])

    workflow.snakemake_cmd_absolute = harpy_config["snakemake"]["absolute"]
    workflow.snakemake_cmd_relative = harpy_config["snakemake"]["relative"]
    workflow.config = harpy_config
    workflow.start_text = workflow_info(
        ("Workflow:", workflow.name.replace("_", " ")),
        ("Output Folder:", directory + "/")
    )

    workflow.write_workflow_config()

    workflow.print_onstart()
    workflow.launch(absolute)
#