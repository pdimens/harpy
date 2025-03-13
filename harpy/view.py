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


@click.command(context_settings=dict(allow_interspersed_args=False))
@click.option('-s', '--snakefile',  is_flag = True, show_default = True, default = False, help = "View the snakefile instead")
@click.option('-c', '--config',  is_flag = True, show_default = True, default = False, help = "View the workflow config file instead")
@click.argument('directory', required=True, type=click.Path(exists=True, file_okay=False))
def view(directory, snakefile, config):
    """
    View a workflow log, config, or snakefile

    This convenience command lets you view the latest workflow log file
    of a Harpy output directory. Use `--snakefile` or `--config` to view the workflow
    snakefile or config.yaml file instead, respectively. Use the typical `less` keyboard
    bindings to navigate the output, e.g.:
    
    | key                     | function                   |
    | :---------------------- | :------------------------- |
    | `Up/Down arrow`         | scroll up/down             |
    | `Page Up/Down`          | scroll up/down, but faster |
    | `/` + `pattern`         | search for `pattern`       |
    | `q`                     | exit                       |
    """
    # check if there is a workflow or log folder
    # and whether the expected files are in there
    if snakefile and config:
        print_error("invalid options", "Please pick one of [bold]--snakefile[/bold] or [bold]--config[/bold]")
        sys.exit(1)
    err = 0
    if snakefile:
        files = [i for i in glob.iglob(f"{directory}/workflow/*.smk")]
        err_dir = f"{directory}/workflow/"
        err_file = "There are no snakefiles"
        if not os.path.exists(f"{directory}/workflow"):
            err = 1
        elif not files:
            err = 2
    elif config:
        files = [f"{directory}/workflow/config.yaml"]
        err_dir = f"{directory}/workflow/"
        err_file = "There is no [blue]config.yaml[/blue] file"
        if not os.path.exists(f"{directory}/workflow"):
            err = 1
        elif not os.path.exists(f"{directory}/workflow/config.yaml"):
            err = 2
    else:
        files = [i for i in glob.iglob(f"{directory}/logs/snakemake/*.log*")]
        err_dir = f"{directory}/logs/snakemake/"
        err_file = "There are no log files"
        if not os.path.exists(f"{directory}/logs/snakemake"):
            err = 1
        elif not files:
            err = 2
    if err == 1:
        print_error(
            "directory not found", 
            f"The file you are trying to view is expected to be in [blue]{err_dir}[/blue], but that directory was not found. Please check that this is the correct folder."
        )
        sys.exit(1)
    elif err == 2:
        print_error(
            "file not found", 
            f"{err_file} in [blue]{err_dir}[/blue]. Please check that this is the correct folder."
        )
        sys.exit(1)
    # sort and pull only the most recent file (based on modification time)
    file = sorted(files, key = os.path.getmtime)[-1]
    if not os.access(file, os.R_OK):
        print_error(
            "incorrect permissions",
            f"[blue]{file}[/blue] does not have read access. Please check the file permissions."
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
    click.echo_via_pager(_read_file(file), color = n_colors > 0)
    rprint(
        Panel(
            file,
            title = "[bold blue] File viewed",
            title_align = "left",
            border_style = "dim",
            width = 75
            ),
        file = sys.stderr
    )