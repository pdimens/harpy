"""View the latest log, config, or snakefile of a workflow"""

import os
import sys
import glob
import gzip
import curses
from pygments import highlight
from pygments.lexers import YamlLexer
from click import echo_via_pager
import rich_click as click
from rich.panel import Panel
from rich import print as rprint
from ._printing import print_error
from ._validations import is_gzip

def check_terminal_colors():
    # Initialize curses
    stdscr = curses.initscr()
    # Check if the terminal supports colors
    if not curses.has_colors():
        curses.endwin()
        return 0
    # Start color mode
    curses.start_color()
    # Get the number of colors supported
    num_colors = curses.COLORS
    # Determine the color type based on the number of colors
    if num_colors <= 8:
        ncol = 8
    elif num_colors == 256:
        ncol = 256
    elif num_colors > 256:  # Truecolor (24-bit)
        ncol = 999
    else:
        ncol = 256
    curses.endwin()
    return ncol

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
        sys.exit(1)
    n_colors = check_terminal_colors()
    if n_colors <= 8:
        from pygments.formatters import TerminalFormatter
        formatter = TerminalFormatter
    elif n_colors == 256:
        from pygments.formatters import Terminal256Formatter
        formatter = Terminal256Formatter
    else:
        from pygments.formatters import TerminalTrueColorFormatter
        formatter = TerminalTrueColorFormatter

    def _read_file(x):
        compressed = is_gzip(x)
        opener = gzip.open if compressed else open
        mode = "rt" if compressed else "r"
        with opener(x, mode) as f:
            for line in f:
                yield highlight(line, YamlLexer(),formatter())
    os.environ["PAGER"] = "less -R"
    click.echo_via_pager(_read_file(infile), color = n_colors > 0)
    return infile

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def view():
    """
    View a workflow's components

    These convenient commands let you view the latest workflow log file, snakefile, snakemake parameter
    file, or workflow config file in a directory that was used for the output of a Harpy run.
    Use the typical `less` keyboard bindings to navigate the output, e.g.:
    
    | key                     | function                   |
    | :---------------------- | :------------------------- |
    | `Up/Down` arrow         | scroll up/down             |
    | `Page Up/Down`          | faster up/down scrolling   |
    | `/` + `pattern`         | search for `pattern`       |
    | `q`                     | exit                       |
    """

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False))
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False), nargs=1)
def config(directory):
    """
    View a workflow's config file
    
    The workflow config file has all of the parameters and user inputs that went into the workflow.
    The only required input is the output folder designated in a previous Harpy run, where you can find
    `workflow/config.harpy.yaml`. Navigate with the typical `less` keyboard bindings, e.g.:
    
    | key                     | function                   |
    | :---------------------- | :------------------------- |
    | `Up/Down` arrow         | scroll up/down             |
    | `Page Up/Down`          | faster up/down scrolling   |
    | `/` + `pattern`         | search for `pattern`       |
    | `q`                     | exit                       |
    """
    target_file = f"{directory}/workflow/config.harpy.yaml"
    err_dir = f"{directory}/workflow/"
    err_file = "There is no [blue]config.harpy.yaml[/] file"
    if not os.path.exists(f"{directory}/workflow"):
        print_error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/], but that directory was not found. Please check that this is the correct folder."
        )
        sys.exit(1)
    elif not os.path.exists(target_file):
        print_error(
            "file not found", 
            f"{err_file} in [blue]{err_dir}[/]. Please check that this is the correct folder."
        )
        sys.exit(1)
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

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False))
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False), nargs=1)
def log(directory):
    """
    View a workflow's last log file
    
    The log file contains everything Snakemake printed during runtime.
    The only required input is the output folder designated in a previous Harpy run, where you can find
    `logs/snakemake/`. Navigate with the typical `less` keyboard bindings, e.g.:
    
    | key                     | function                   |
    | :---------------------- | :------------------------- |
    | `Up/Down` arrow         | scroll up/down             |
    | `Page Up/Down`          | faster up/down scrolling   |
    | `/` + `pattern`         | search for `pattern`       |
    | `q`                     | exit                       |
    """
    files = [i for i in glob.iglob(f"{directory}/logs/snakemake/*.log*")]
    target_file = sorted(files, key = os.path.getmtime)[-1]
    err_dir = f"{directory}/logs/snakemake/"
    err_file = "There are no log files"
    if not os.path.exists(f"{directory}/logs/snakemake"):
        print_error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/], but that directory was not found. Please check that this is the correct folder."
        )
        sys.exit(1)
    elif not files:
        print_error(
            "files not found", 
            f"{err_file} in [blue]{err_dir}[/]. Please check that this is the correct folder."
        )
        sys.exit(1)
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

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False))
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False), nargs=1)
def snakefile(directory):
    """
    View a workflow's snakefile
    
    The snakefile contains all the instructions for a workflow.
    The only required input is the output folder designated in a previous Harpy run, where you can find
    `workflow/<workflow>.smk`. Navigate with the typical `less` keyboard bindings, e.g.:
    
    | key                     | function                   |
    | :---------------------- | :------------------------- |
    | `Up/Down` arrow         | scroll up/down             |
    | `Page Up/Down`          | faster up/down scrolling   |
    | `/` + `pattern`         | search for `pattern`       |
    | `q`                     | exit                       |
    """
    files = [i for i in glob.iglob(f"{directory}/workflow/*.smk")]
    target_file = files[0]
    err_dir = f"{directory}/workflow/"
    err_file = "There are no snakefiles"
    if not os.path.exists(f"{directory}/workflow"):
        print_error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/], but that directory was not found. Please check that this is the correct folder."
        )
        sys.exit(1)
    elif not files:
        print_error(
            "file not found", 
            f"{err_file} in [blue]{err_dir}[/]. Please check that this is the correct folder."
        )
        sys.exit(1)
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

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False))
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False), nargs=1)
def snakeparams(directory):
    """
    View a workflow's snakemake parameter file
    
    The snakemake parameter file file has the runtime parameters snakemake was invoked with.
    The only required input is the output folder designated in a previous Harpy run, where you can find
    `workflow/config.yaml`. Navigate with the typical `less` keyboard bindings, e.g.:
    
    | key                     | function                   |
    | :---------------------- | :------------------------- |
    | `Up/Down` arrow         | scroll up/down             |
    | `Page Up/Down`          | faster up/down scrolling   |
    | `/` + `pattern`         | search for `pattern`       |
    | `q`                     | exit                       |
    """
    target_file = f"{directory}/workflow/config.yaml"
    err_dir = f"{directory}/workflow/"
    err_file = "There is no [blue]config.yaml[/] file"
    if not os.path.exists(f"{directory}/workflow"):
        print_error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/], but that directory was not found. Please check that this is the correct folder."
        )
        sys.exit(1)
    elif not os.path.exists(target_file):
        print_error(
            "file not found", 
            f"{err_file} in [blue]{err_dir}[/]. Please check that this is the correct folder."
        )
        sys.exit(1)
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
view.add_command(log)
view.add_command(snakefile)
view.add_command(snakeparams)