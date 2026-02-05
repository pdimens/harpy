"""Module to bypass Harpy and run snakemake"""

from datetime import datetime
import os
import sys
import yaml
import rich_click as click
from harpy.common.environments import check_environments
from harpy.common.printing import print_error, workflow_info
from harpy.common.workflow import Workflow

def config_extract(d: dict, section:str, key: str = "", allow_missing: bool = False):
    val = d.get(section, {})
    if not val:
        if not allow_missing:
            print_error(
                "incorrect workflow.yaml",
                f"The [blue]workflow.yaml[/] file is missing the [green]{section}[/] section",
                f"Please verify that [blue]workflow/workflow.yaml[/] is not missing the [green]{section}[/] section (and it is not empty)."
            )
        else:
            return val
    if ":" in key:
        _key, _subkey = key.split(":")
        val = val.get(_key, {})
        if val:
            val = val.get(_subkey, {})
    elif key:
        val = val.get(key, {})

    if not val:
        print_error(
            "incorrect workflow.yaml",
            "The [blue]workflow.yaml[/] file is missing one or more the necessary/expected keys.",
            f"Please verify that [blue]workflow/workflow.yaml[/] is not missing the [green]{key}[/] key (and it is not empty) under the [green]{section}[/] section."
        )
    return val

def snakemake_profile_extract(d: dict, key: str):
    val = d.get(key, None)
    if not val:
        print_error(
            "incorrect config.yaml",
            "The [blue]config.yaml[/] file is missing one or more the necessary/expected keys.",
            f"Please verify that [blue]workflow/config.yaml[/] is not missing the [green]{key}[/] key (and it is not empty)."
        )
    return val

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/other")
@click.option('-a', '--absolute',  is_flag = True, default = False, help = 'Call Snakemake with absolute paths')
@click.option('-d', '--direct',  is_flag = True, default = False, help = 'Call Snakemake directly without Harpy intervention')
@click.option('-t', '--threads', type = click.IntRange(2, 999, clamp = True), help = 'Change the number of threads (>1)')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--quiet', default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.help_option('--help', panel = "Workflow Options", hidden = True)
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
    _name = config_extract(harpy_config, "Workflow", "name")
    _allow_noparams = True if "validate" in _name else False
    _dir = snakemake_profile_extract(snakemake_config, "directory")
    _inputs = config_extract(harpy_config, "Inputs")

    workflow = Workflow(_name, "NA", _dir, False, clean, quiet)
    if isinstance(_inputs, list):
        workflow.input(_inputs)
    else:
        for i,j in _inputs.items():
            workflow.input(j, i)    

    workflow.parameters = config_extract(harpy_config, "Parameters", allow_missing=_allow_noparams)

    workflow.snakemake_cmd_relative = config_extract(harpy_config, "Workflow", "snakemake:relative")
    if absolute:
        workflow.snakemake_cmd_absolute = config_extract(harpy_config, "Workflow", "snakemake:absolute")
    else:
        _abs_cmd = []
        for i in workflow.snakemake_cmd_relative.split():
            if os.path.exists(i):
                _abs_cmd.append(os.path.abspath(i))
            else:
                _abs_cmd.append(i)
        workflow.snakemake_cmd_absolute = " ".join(_abs_cmd)

    _sdm = snakemake_profile_extract(snakemake_config, "software-deployment-method")
    if _sdm != "apptainer":
        workflow.conda = config_extract(harpy_config, "Workflow", "snakemake:conda-envs")
        check_environments(_dir, workflow.conda)

    # inherit workflow report part, if present
    workflow.notebooks = harpy_config["Workflow"].get("reports", {"skip" : False})

    # inherit workflow report part, if present
    workflow.linkedreads = harpy_config["Workflow"].get("linkedreads", {"type" : "none"})

    _smlog = os.path.join(_dir, config_extract(harpy_config, "Workflow", "snakemake:log"))

    # overwrite workflow.yaml with the new logfile name, otherwise just create the
    # workflow.config dict
    if os.path.exists(_smlog) or os.path.exists(_smlog + ".gz"):
        timestamp = datetime.now().strftime("%d_%m_%Y") + ".log"
        split_log = _smlog.split(".")
        _basename = ".".join(split_log[0:-3])
        incremenent = int(split_log[-3]) + 1
        workflow.snakemake_logfile = f"{_basename}.{incremenent}.{timestamp}"

    workflow.write_workflow_config()

    if threads:
        snakemake_config["cores"] = threads
        workflow.profile = snakemake_config
        workflow.write_snakemake_profile()

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
