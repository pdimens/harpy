"""Module to bypass Harpy and run snakemake"""

from ast import keyword
from datetime import datetime
import os
import re
import sys
import yaml
import rich_click as click
from harpy.common.conda import check_environments
from harpy.common.printing import print_error, workflow_info
from harpy.common.workflow import Workflow

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/other")
@click.option('-a', '--absolute',  is_flag = True, default = False, help = 'Call Snakemake with absolute paths')
@click.option('-d', '--direct',  is_flag = True, default = False, help = 'Call Snakemake directly without Harpy intervention')
@click.option('-t', '--threads', type = click.IntRange(2, 999, clamp = True), help = 'Change the number of threads (>1)')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--quiet', default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False, readable=True, resolve_path=True), nargs=1)
def resume(directory, absolute, direct, threads, clean, quiet):
    """
    Continue an incomplete Harpy workflow

    In the event you need to run the Snakemake workflow present in a Harpy output directory
    (e.g. `Align/bwa`) without redoing validations and writing new configuration files (config files get updated with new thread count and log file name),
    this command bypasses all the preprocessing steps of Harpy workflows and executes the Snakemake command
    present in `directory/workflow/workflow.yaml`.

    The target directory must have:
    - the `workflow/config.yaml` file
    - the `workflow/workflow.yaml` file
    - `workflow/envs/*.yaml` file(s) if using conda
    - `workflow/hpc/config.yaml` if using HPC
    """
    CONFIG_FILE = os.path.join(directory, "workflow", "workflow.yaml")
    PROFILE_FILE = os.path.join(directory, "workflow", "config.yaml")
    if not os.path.exists(PROFILE_FILE):
        print_error("missing snakemake config", f"Target directory [yellow]{directory}[/] does not contain the file [blue]workflow/config.yaml[/]")
    if not os.path.exists(CONFIG_FILE):
        print_error("missing workflow config", f"Target directory [yellow]{directory}[/] does not contain the file [blue]workflow/workflow.yaml[/]")
    
    with open(CONFIG_FILE, 'r', encoding="utf-8") as f:
        harpy_config: dict = yaml.full_load(f)
    with open(PROFILE_FILE, 'r', encoding="utf-8") as f:
        snakemake_config: dict = yaml.full_load(f)

    #container = snakemake_config["software-deployment-method"] == "apptainer"
    try:
        workflow = Workflow(harpy_config["workflow"], "NA", snakemake_config["directory"], False, clean, quiet)
        # pull in the inputs and store them, removing the original
        workflow.inputs = harpy_config.pop("inputs")
    except KeyError:
        print_error(
            "incorrect config.yaml",
            "The [blue]workflow.yaml[/] file is missing one or more the necessary/expected keys.",
            f"Please verify that [blue]{os.path.relpath(directory)}/workflow/workflow.yaml[/] is not missing these keys (and they are not empty): [green]workflow[/], [green]snakemake[/], [green]inputs[/]"
        )

    workflow.config = harpy_config
    
    if snakemake_config["software-deployment-method"] != "apptainer":
        try:
            workflow.conda = harpy_config["snakemake"]["conda_envs"]
        except KeyError:
            print_error(
                "missing conda environments",
                "The [blue]config.yaml[/] file indicates Snakemake will run with conda environments (rather than apptainer), but the [blue]workflow.yaml[/] file is missing the [green]snakemake:conda-envs[/] key.",
                f"If you intend to run Harpy with Snakemake using conda environments, [blue]{os.path.relpath(directory)}/workflow/workflow.yaml[/] will need [green]snakemake:conda-envs[/] specified. That file can be rebuilt using [green]harpy {workflow.name} --setup ..."
            )   
        check_environments(directory, harpy_config["snakemake"]["conda_envs"])

    try:
        sm_log = os.path.join(directory, harpy_config["snakemake"]["log"])
    except KeyError:
        print_error(
            "incorrect config.yaml",
            "The [blue]workflow.yaml[/] file is missing the [green]snakemake:log[/] key.",
            f"Please verify that [blue]{os.path.relpath(directory)}/workflow/workflow.yaml[/] is not missing this key and it's not empty."
        )

    if os.path.exists(sm_log) or os.path.exists(sm_log + ".gz"):
        timestamp = datetime.now().strftime("%d_%m_%Y") + ".log"
        split_log = sm_log.split(".")
        _basename = ".".join(split_log[0:-3])
        incremenent = int(split_log[-3]) + 1
        workflow.snakemake_logfile = f"{_basename}.{incremenent}.{timestamp}"
        harpy_config["snakemake"]["log"] = workflow.snakemake_logfile
        workflow.write_workflow_config()

    if threads:
        snakemake_config["cores"] = threads
        workflow.profile = snakemake_config
        workflow.write_snakemake_profile()

    try:
        sm_abs = harpy_config["snakemake"].get("absolute", None)
        workflow.snakemake_cmd_absolute = sm_abs
        if not sm_abs and absolute:
            print_error(
                "missing Snakemake command",
                "You requested the Snakemake absolute-path command but it was missing from [blue]workflow.yaml[/]. Please check that the configuration has the [green]snakemake:absolute[/] key and is associated with a Snakemake command line call."
            )
    except KeyError:
        print_error(
            "incorrect config.yaml",
            "The [blue]workflow.yaml[/] file is missing the [green]snakemake[/] key.",
            f"Please verify that [blue]{os.path.relpath(directory)}/workflow/workflow.yaml[/] is not missing this key and it's not empty."
        )
    
    try:
        sm_rel = harpy_config["snakemake"].get("relative", None)
        workflow.snakemake_cmd_relative = sm_rel
        if not sm_rel and not absolute:
            print_error(
                "missing Snakemake command",
                "Attempted to use the Snakemake relative-path command but it was missing from [blue]workflow.yaml[/]. Please check that the configuration has the [green]snakemake:relative[/] key and is associated with a Snakemake command line call."
            )
    except KeyError:
        print_error(
            "incorrect config.yaml",
            "The [blue]workflow.yaml[/] file is missing the [green]snakemake[/] key.",
            f"Please verify that [blue]{os.path.relpath(directory)}/workflow/workflow.yaml[/] is not missing this key and it's not empty."
        )

    workflow.start_text = workflow_info(
        ("Workflow:", workflow.name.replace("_", " ")),
        ("Output Folder:", directory + "/")
    )

    if direct:
        if absolute:
            _ = os.system(workflow.snakemake_cmd_absolute)
        else:
            _ = os.system(workflow.snakemake_cmd_relative)
        if _ > 0:
            sys.exit(1)
    else:
        workflow.print_onstart()
        workflow.launch(absolute)
