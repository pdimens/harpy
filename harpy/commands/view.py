"""View the latest log, config, or snakefile of a workflow"""

import os
import glob
import rich_click as click
from rich.panel import Panel
from rich.tree import Tree
from harpy.common.file_ops import choose_logfile, parse_error, parse_file
from harpy.common.printing import HarpyPrint

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
    Browse or edit a workflow's config file
    
    The workflow config file has all of the parameters and user inputs that went into the workflow.
    The only required input is the output folder designated in a previous Harpy run, where you can find
    `workflow/workflow.yaml`.
    """
    hp = HarpyPrint()
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
    hp.print(
        Panel(
            target_file,
            title = "[bold blue] File viewed",
            title_align = "left",
            border_style = "dim",
            width = 75
        )
    )

@click.command()
@click.help_option('--help', hidden = True)
@click.argument('program', required=False, type=str)
def envs(program):
    """
    Print the Snakemake-managed conda environments

    This convenience command will print the main information of the conda environment recipes within
    `.environments/`, which can be useful when troubleshooting requires you to enter a specific conda environment.
    Optionally provide the name of a `PROGRAM` (or partial name, case insensitive) to only return the environments that
    contain that program (e.g. `leviath` would match `leviathan`).
    """
    hp = HarpyPrint()
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
    tree = Tree("[bold]Conda Environments")
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
        stripdeps = [j.split("=")[0] for j in deps]
        if not program or any([program in j for j in stripdeps]):
            _subtree = tree.add("[bold]" + i.removesuffix('.yaml'), style = "blue")
            for d in deps:
                if program:
                    if program.lower() in d:
                        _subtree.add(d + " <", style = 'bold green', highlight = False)
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
    Browse a workflow's Snakemake log
    
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
    hp = HarpyPrint()
    target_file = choose_logfile(directory, choose)
    parse_file(target_file)
    hp.print(
        Panel(
            target_file,
            title = "[bold blue] File viewed",
            title_align = "left",
            border_style = "dim",
            width = 75
        )
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
    Browse or edit a workflow's Snakefile
    
    The snakefile contains all the instructions for a workflow. The only required input is the output folder
    previously created by Harpy where you can find `workflow/workflow.smk`.
    """
    hp = HarpyPrint()
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
    hp.print(
        Panel(
            target_file,
            title = "[bold blue] File viewed",
            title_align = "left",
            border_style = "dim",
            width = 75
        )
    )

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option("-e", "--edit", is_flag=True, default=False, help = "Open the config file in you system's default editor")
@click.help_option('--help', hidden = True)
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False), nargs=1)
def profile(directory, edit):
    """
    Browse or edit a workflow's Snakemake configurations
    
    The snakemake configuration file has the runtime parameters snakemake was invoked with (i.e.,
    computational specifics that don't impact your results). The only required input is the output folder
    previously created by Harpy where you can find `workflow/profile.yaml`.
    """
    hp = HarpyPrint()
    err_dir = os.path.join(directory, "workflow")
    target_file = os.path.join(err_dir, "profile.yaml")
    err_file = "There is no [blue]profile.yaml[/] file"
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
    hp.print(
        Panel(
            target_file,
            title = "[bold blue] File viewed",
            title_align = "left",
            border_style = "dim",
            width = 75
        )
    )

view.add_command(config)
view.add_command(envs)
view.add_command(error)
view.add_command(log)
view.add_command(snakefile)
view.add_command(profile)