"""View the latest log, config, or snakefile of a workflow"""

import os
import sys
import glob
import rich_click as click
from rich.panel import Panel
from rich import print as rprint
from rich.tree import Tree
from harpy.common.file_ops import choose_logfile, parse_error, parse_file
from harpy.common.printing import HarpyPrint

hp = HarpyPrint()

@click.group(options_metavar='')
@click.help_option('--help', hidden = True)
def view():
    """
    View a workflow's components

    These convenient commands let you view/edit the latest workflow log file, snakefile, snakemake parameter
    file, workflow config file in a directory that was used for the output of a Harpy run.
    """

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option("-e", "--edit", is_flag=True, default=False, help = "Open the config file in you system's default editor")
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
@click.help_option('--help', hidden = True)
def config(directory, edit):
    """
    View/edit a workflow's config file
    
    The workflow config file has all of the parameters and user inputs that went into the workflow.
    The only required input is the output folder designated in a previous Harpy run, where you can find
    `workflow/workflow.yaml`.
    """
    err_dir = os.path.join(directory, "workflow")
    target_file = os.path.join(err_dir, "workflow.yaml")
    err_file = "There is no [blue]workflow.yaml[/] file"
    if not os.path.exists(err_dir):
        hp.error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/], but that directory was not found. Please check that this is the correct folder."
        )
    elif not os.path.exists(target_file):
        hp.error(
            "file not found", 
            f"{err_file} in [blue]{err_dir}[/]. Please check that this is the correct folder."
        )
    if edit:
        click.edit(filename = target_file, extension = "yaml")
    else:
        parse_file(target_file)
    rprint(
        Panel(
            target_file,
            title = "[bold blue] File viewed",
            title_align = "left",
            border_style = "dim",
            width = 75
            ),
        file = sys.stderr
    )

@click.command()
@click.help_option('--help', hidden = True)
@click.argument('program', required=False, type=str)
def environments(program):
    """
    Print the Snakemake-managed conda environments

    This convenience command will print the main information of the conda environment recipes within
    `.environments/`, which can be useful when troubleshooting requires you to enter a specific conda environment.
    Optionally provide the name of a `PROGRAM` (or partial name, case insensitive) to only return the environments that
    contain that program (e.g. `leviath` would match `leviathan`).
    """
    if not os.path.exists(".environments"):
        hp.error(
            "directory not found", 
            "No [blue].environments/[/] folder found in the current directory."
        )
    files = [i for i in glob.iglob(".environments/*.yaml")]
    if not files:
        hp.error(
            "files not found", 
            "No conda recipes ending in [green].yaml[/] found in [blue].environments[/]."
        )
    tree = Tree("[bold light_steel_blue]Conda Environments")
    for i in files:
        deps = []
        with open(i, "r") as file:
            skip = True
            for line in file:
                if line.startswith("dependencies"):
                    skip = False
                    continue
                if not skip:
                    dep = line.split("::")[-1].rstrip()
                    deps.append(dep.rstrip())
                    #deps += f" {dep.rstrip()}"
        if (program and program.lower() in deps) or not program:
            _subtree = tree.add(i.removesuffix('.yaml'), style = "bold blue")
            for d in deps:
                if program:
                    if program.lower() in d:
                        _subtree.add(d, style = 'bold blue', highlight = False)
                    else:
                        _subtree.add(d, style = 'dim default', highlight = False)
                else:
                    _subtree.add(d, style = 'default', highlight = False)

    hp.print(tree)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option("-c", "--choose", is_flag=True, default=False, help = "List logs for user choice")
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
@click.help_option('--help', hidden = True)
def log(directory, choose):
    """
    View a workflow's Snakemake log file
    
    The log file contains everything Snakemake printed during runtime.
    The only required input is an output folder created by Harpy where you can find
    `.snakemake/log`. Use `--choose` to pick from a list of all Snakemake logfiles in
    the `directory`. Navigate with the typical `less` keyboard bindings, e.g.:
    
    | key                     | function                   |
    | :---------------------- | :------------------------- |
    | `Up/Down` arrow         | scroll up/down             |
    | `Page Up/Down`          | faster up/down scrolling   |
    | `/` + `pattern`         | search for `pattern`       |
    | `q`                     | exit                       |
    """
    target_file = choose_logfile(directory, choose)
    parse_file(target_file)
    rprint(
        Panel(
            target_file,
            title = "[bold blue] File viewed",
            title_align = "left",
            border_style = "dim",
            width = 75
            ),
        file = sys.stderr
    )

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option("-c", "--choose", is_flag=True, default=False, help = "List logs for user choice")
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
@click.help_option('--help', hidden = True)
def error(directory, choose):
    """
    Print a workflow's error
    
    The log file contains everything Snakemake printed during runtime and this command
    scans the file for an error and prints it to the terminal. The only required input
    is an output folder created by Harpy where you can find `.snakemake/log`. Use 
    `--choose` to pick from a list of all Snakemake logfiles in the `directory`.
    """
    target_file = choose_logfile(directory, choose)
    parse_error(target_file)


@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option("-e", "--edit", is_flag=True, default=False, help = "Open the config file in you system's default editor")
@click.help_option('--help', hidden = True)
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
def snakefile(directory, edit):
    """
    View/edit a workflow's Snakefile
    
    The snakefile contains all the instructions for a workflow. The only required input is the output folder
    previously created by Harpy where you can find `workflow/workflow.smk`.
    """
    workdir = os.path.join(directory, "workflow")
    if not os.path.exists(workdir):
        hp.error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{workdir}[/], but that directory was not found. Please check that you are looking in the correct folder."
        )

    target_file = os.path.join(workdir, "workflow.smk")
    if not os.path.exists(target_file):
        hp.error(
            "snakefile not found", 
            f"[blue]{target_file}[/] was not found. Please check that you are looking in the correct folder."
        )
    if edit:
        click.edit(filename = target_file, extension = "yaml")
    else:
        parse_file(target_file)
    rprint(
        Panel(
            target_file,
            title = "[bold blue] File viewed",
            title_align = "left",
            border_style = "dim",
            width = 75
            ),
        file = sys.stderr
    )

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option("-e", "--edit", is_flag=True, default=False, help = "Open the config file in you system's default editor")
@click.help_option('--help', hidden = True)
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False), nargs=1)
def snakeparams(directory, edit):
    """
    View/edit a workflow's Snakemake configurations
    
    The snakemake configuration file has the runtime parameters snakemake was invoked with (i.e.,
    computational specifics that don't impact your results). The only required input is the output folder
    previously created by Harpy where you can find `workflow/config.yaml`.
    """
    err_dir = os.path.join(directory, "workflow")
    target_file = os.path.join(err_dir, "config.yaml")
    err_file = "There is no [blue]config.yaml[/] file"
    if not os.path.exists(err_dir):
        hp.error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/], but that directory was not found. Please check that you are looking in the correct folder."
        )
    elif not os.path.exists(target_file):
        hp.error(
            "file not found", 
            f"{err_file} in [blue]{err_dir}[/]. Please check that you are looking in the the correct folder."
        )
    if edit:
        click.edit(filename = target_file, extension = "yaml")
    else:
        parse_file(target_file)
    rprint(
        Panel(
            target_file,
            title = "[bold blue] File viewed",
            title_align = "left",
            border_style = "dim",
            width = 75
            ),
        file = sys.stderr
    )

view.add_command(config)
view.add_command(environments)
view.add_command(error)
view.add_command(log)
view.add_command(snakefile)
view.add_command(snakeparams)