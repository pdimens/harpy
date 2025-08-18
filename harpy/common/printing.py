"""Module of pretty-printing for errors and prompts"""

import time as _time
import os
import sys
from rich.console import Console, RenderableType
from rich import box
from rich.table import Table
from rich.panel import Panel

CONSOLE = Console(stderr=True)

def print_error(
    errortitle: str,
    errortext: RenderableType|str,
    solutiontext: RenderableType|str|None = None,
    offendertitle: str|None = None,
    offenders: RenderableType|list|str|None = None,
    _exit: bool = True
    ) -> None:
    """
    Print a yellow panel with error text, exits with error code 1 by default, but exit
    can be disabled with _exit = False. Prints a blue-panel solution if `solutiontext` is provided,
    prints `offenders` after solution if proivded.    
    """
    CONSOLE.print(
        Panel(
            errortext,
            title = f"[bold]Error: {errortitle}",
            title_align = "left",
            border_style = "yellow",
            width = 75
            )
    )
    if solutiontext:
        CONSOLE.print(
            Panel(solutiontext,
                title = "[bold]Solution",
                title_align = "left",
                border_style = "blue",
                subtitle= offendertitle,
                width = 75
            )
        )
        if offenders:
            offenders = offenders if isinstance(offenders, str) else "\n".join([str(i) for i in offenders])
            CONSOLE.print(
                offenders,
                sep = "\n",
                highlight= False
            )
    if _exit:
        sys.exit(1)

def print_notice(noticetext: str) -> None:
    """Print a white panel with information text text"""
    CONSOLE.print(
        Panel(
            noticetext,
            title = "[dim]Notice",
            title_align = "left",
            border_style = "dim",
            width = 75
            )
    )

def print_onstart(text: str, title: str) -> None:
    """Print a panel of info on workflow run"""
    CONSOLE.print("")
    CONSOLE.rule(f"[bold]harpy {title}", style = "light_steel_blue")
    CONSOLE.print(text)

def print_setup_error(exitcode: int) -> None:
    """Print a red panel with snakefile or conda/singularity error text"""
    if exitcode == 1:
        errortype = "Snakefile Error"
        errortext = "Something is wrong with the Snakefile for this workflow. If you manually edited the Snakefile, see the error below for troubleshooting. If you didn't, it's probably a bug (oops!) and you should submit an issue on GitHub: [bold]https://github.com/pdimens/harpy/issues"
    else:
        errortype = "Software Environment Error"
        errortext = "There was an issue creating the software environment necessary to run this workflow. If you manually edited the conda dependencies in [blue]/workflows/envs[/], see the error below for troubleshooting. If you didn't, it might be a bug or related to how your system is setup for Conda or Singularity environments and you should submit an issue on GitHub: [bold]https://github.com/pdimens/harpy/issues"
        # Check if this is the `base` conda environment
        current_env = os.environ.get('CONDA_DEFAULT_ENV')
        if current_env == 'base':
            errortext += "\n[yellow]Notice:[/] Harpy detected that you're in the [blue]base[/] conda environment; conda recommends against installing anything into the [base]base[/] environment."
        # Check the channel priority setting
        try:
            from conda.base.context import context
            if context.channel_priority == "strict":
                errortext += "\n[yellow]Notice:[/] Your conda channel priority is configured as [yellow]strict[/], which can sometimes cause issues with Snakemake creating conda environments. Ignore this detail if you are using [blue]--container[/]."
        except ModuleNotFoundError:
            pass
    CONSOLE.rule(f"[bold]{errortype}[/] [dim]" + _time.strftime('%d %b %Y @ %H:%M'), style = "red")
    CONSOLE.print(errortext)
    CONSOLE.rule("[bold]Error Reported by Snakemake", style = "red")

def print_onerror(logfile: str, time = None) -> None:
    """Print a red panel with error text. To be used in place of onerror: inside a snakefile. Expects the erroring rule printed after it."""
    days = time.days
    seconds = time.seconds
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    seconds = seconds % 60
    time_text = f"{days} days, {hours} hours, {minutes} minutes, {seconds} seconds"
    datatable = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    datatable.add_column("detail", justify="left", style="red", no_wrap=True)
    datatable.add_column("value", justify="left")
    datatable.add_row("Duration:", time_text)
    datatable.add_row("Workflow Log: ", logfile + ".gz")
    CONSOLE.rule("[bold]Workflow Error[/][dim] " + _time.strftime('%d %b %Y @ %H:%M'), style = "red")
    CONSOLE.print("The workflow stopped because of an error. See the information Snakemake reported below.")
    CONSOLE.print(datatable)
    CONSOLE.rule("[bold]Where Error Occurred", style = "red")

def workflow_info(*arg: tuple[str, str | int | float]) -> Table:
    """
    Accepts an unlimited number of length-2 lists or tuples and returns a rich.Table with the value of the first indices as the row names and the second indices as the values
    Use None instead of a list to ignore that entry (useful for conditionals). The second value will always be converted to a string.
    """
    table = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    table.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    table.add_column("value", justify="left")
    for i in arg:
        if i:
            table.add_row(i[0], str(i[1]))
    return table