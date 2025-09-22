"""View the latest log, config, or snakefile of a workflow"""

import os
import sys
import glob
import gzip
import curses
from pygments import highlight
from pygments.lexers import get_lexer_by_name
from pygments.formatters import get_formatter_by_name
from click import echo_via_pager
import rich_click as click
from rich.console import Console
from rich.panel import Panel
from rich import print as rprint
from harpy.common.printing import print_error
from harpy.common.file_ops import is_gzip

def check_terminal_colors():
    # Initialize curses and always tear down
    try:
        _ = curses.initscr()
        # Check if the terminal supports colors
        if not curses.has_colors():
            return 0
        curses.start_color()
        num_colors = curses.COLORS
        # Determine the color type based on the number of colors
        if num_colors <= 8:
            return 8
        else:
            return 256
    except curses.error:
        # Non-interactive/unsupported terminals
        return 0
    finally:
        try:
            curses.endwin()
        except curses.error:
            pass

def parse_file(infile):
    '''
    take a list of input file name, get the most recent by modificiation time, and print it via pygmentized less
    returns a string of the file that was viewed
    '''
    if not os.access(infile, os.R_OK):
        print_error(
            "incorrect permissions",
            f"[blue]{infile}[/] does not have read access. Please check the file permissions."
        )
    n_colors = check_terminal_colors()
    if n_colors <= 8:
        formatter = get_formatter_by_name("terminal")
    else:
        formatter = get_formatter_by_name("terminal256")

    def _read_file(x):
        compressed = is_gzip(x)
        opener = gzip.open if compressed else open
        mode = "rt" if compressed else "r"
        lexer = get_lexer_by_name("yaml")
        with opener(x, mode) as f:
            for line in f:
                yield highlight(line, lexer, formatter)
    os.environ["PAGER"] = "less -R"
    echo_via_pager(_read_file(infile), color = n_colors > 0)
    return infile

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def view():
    """
    View a workflow's components

    These convenient commands let you view/edit the latest workflow log file, snakefile, snakemake parameter
    file, workflow config file in a directory that was used for the output of a Harpy run.
    """

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.option("-e", "--edit", is_flag=True, default=False, help = "Open the config file in you system's default editor")
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False), nargs=1)
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
        print_error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/], but that directory was not found. Please check that this is the correct folder."
        )
    elif not os.path.exists(target_file):
        print_error(
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
@click.argument('program', required=False, type=str, nargs=1)
def environments(program):
    """
    View the Snakemake-managed conda environments

    This convenience command will print the main information of the conda environment recipes within
    `.environments/`, which can be useful when troubleshooting requires you to enter a specific conda environment.
    Optionally provide the name of a `PROGRAM` (or partial name, case insensitive) to only return the environments that
    contain that program (e.g. `leviath` would match `leviathan`).
    """
    if not os.path.exists(".environments"):
        print_error(
            "directory not found", 
            "No [blue].environments/[/] folder found in the current directory."
        )
    files = [i for i in glob.iglob(".environments/*.yaml")]
    if not files:
        print_error(
            "files not found", 
            "No conda recipes ending in [green].yaml[/] found in [blue].environments[/]."
        )
    console = Console()

    for i in files:
        deps = ""
        with open(i, "r") as file:
            skip = True
            for line in file:
                if line.startswith("dependencies"):
                    skip = False
                    continue
                if not skip:
                    dep = line.split("::")[-1].rstrip()
                    deps += f" {dep.rstrip()}"
        if (program and program.lower() in deps) or not program:
            console.print()
            console.rule(i.removesuffix('.yaml'), style = "blue")
            for d in deps.split():
                if program and program.lower() in d:
                    console.print(f"→ {d}", style = "bold green")
                else:
                    console.print(f"- {d}", style = "default", highlight=False)
            console.print()
    return

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False})
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False), nargs=1)
def log(directory):
    """
    View a workflow's last log file
    
    The log file contains everything Snakemake printed during runtime.
    The only required input is the output folder previously created by Harpy where you can find
    `logs/snakemake/`. Navigate with the typical `less` keyboard bindings, e.g.:
    
    | key                     | function                   |
    | :---------------------- | :------------------------- |
    | `Up/Down` arrow         | scroll up/down             |
    | `Page Up/Down`          | faster up/down scrolling   |
    | `/` + `pattern`         | search for `pattern`       |
    | `q`                     | exit                       |
    """
    err_dir = os.path.join(directory, "logs", "snakemake")
    err_file = "There are no log files"
    if not os.path.exists(err_dir):
        print_error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/], but that directory was not found. Please check that this is the correct folder."
        )
    files = [i for i in glob.iglob(f"{err_dir}/*.log*")]        
    if not files:
        print_error(
            "files not found", 
            f"{err_file} in [blue]{err_dir}[/]. Please check that this is the correct folder."
        )
    target_file = sorted(files, key = os.path.getmtime)[-1]
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
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False), nargs=1)
def snakefile(directory, edit):
    """
    View/edit a workflow's snakefile
    
    The snakefile contains all the instructions for a workflow. The only required input is the output folder
    previously created by Harpy where you can find `workflow/workflow.smk`.
    """
    workdir = os.path.join(directory, "workflow")
    if not os.path.exists(workdir):
        print_error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{workdir}[/], but that directory was not found. Please check that you are looking in the correct folder."
        )

    target_file = os.path.join(workdir, "workflow.smk")
    if not os.path.exists(target_file):
        print_error(
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
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False), nargs=1)
def snakeparams(directory, edit):
    """
    View/edit a workflow's snakemake configurations
    
    The snakemake configuration file has the runtime parameters snakemake was invoked with (i.e.,
    computational specifics that don't impact your results). The only required input is the output folder
    previously created by Harpy where you can find `workflow/config.yaml`.
    """
    err_dir = os.path.join(directory, "workflow")
    target_file = os.path.join(err_dir, "config.yaml")
    err_file = "There is no [blue]config.yaml[/] file"
    if not os.path.exists(err_dir):
        print_error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/], but that directory was not found. Please check that you are looking in the correct folder."
        )
    elif not os.path.exists(target_file):
        print_error(
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
view.add_command(log)
view.add_command(snakefile)
view.add_command(snakeparams)