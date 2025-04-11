"""Module of pretty-printing for errors and prompts"""

import sys
from rich import print as rprint
from rich.console import Console
from rich import box
from rich.table import Table
from rich.panel import Panel

console = Console(stderr=True)

def print_error(errortitle: str, errortext: str) -> None:
    """Print a yellow panel with error text"""
    rprint(
        Panel(
            errortext,
            title = f"[bold]Error: {errortitle}",
            title_align = "left",
            border_style = "yellow",
            width = 75
            ),
        file = sys.stderr
    )

def print_solution(solutiontext: str) -> None:
    """Print a blue panel with solution text"""
    rprint(
        Panel(solutiontext,
            title = "[bold]Solution",
            title_align = "left",
            border_style = "blue",
            width = 75
            ),
        file = sys.stderr
    )

def print_solution_with_culprits(solutiontext: str, culprittext: str) -> None:
    """Print a blue panel with solution text and culprittext as the subtitle to introducethe list of offenders below it."""
    rprint(
        Panel(solutiontext,
            title = "[bold]Solution",
            title_align = "left",
            subtitle = culprittext,
            border_style = "blue",
            width = 75
            ),
        file = sys.stderr
    )

def print_notice(noticetext: str) -> None:
    """Print a white panel with information text text"""
    rprint(
        Panel(
            noticetext,
            title = "[dim]Notice",
            title_align = "left",
            border_style = "dim",
            width = 75
            ),
        file = sys.stderr
    )

def print_onstart(text: str, title: str) -> None:
    """Print a panel of info on workflow run"""
    rprint("")
    console.rule(f"[bold]harpy {title}", style = "light_steel_blue")
    console.print(text)

def print_setup_error(exitcode: int) -> None:
    """Print a red panel with snakefile or conda/singularity error text"""
    if exitcode == 1:
        errortext = "Something is wrong with the Snakefile for this workflow. If you manually edited the Snakefile, see the error below for troubleshooting. If you didn't, it's probably a bug (oops!) and you should submit an issue on GitHub: [bold]https://github.com/pdimens/harpy/issues"
        errortype = "Snakefile Error"
    else:
        errortext = "There was an issue creating the software environment necessary to run this workflow. If you manually edited the conda dependencies in [blue]/workflows/envs[/], see the error below for troubleshooting. If you didn't, it might be a bug or related to how your system is setup for Conda or Singularity environments and you should submit an issue on GitHub: [bold]https://github.com/pdimens/harpy/issues"
        errortype = "Software Environment Error"
    console.rule(f"[bold]{errortype}", style = "red")
    console.print(errortext)
    console.rule("[bold]Error Reported by Snakemake", style = "red")

def print_onsuccess(outdir: str, summary = None, time = None) -> None:
    """Print a green panel with success text. To be used in place of onsuccess: inside a snakefile"""
    days = time.days
    seconds = time.seconds
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    seconds = seconds % 60
    time_text = f"{days} days, {hours} hours, {minutes} minutes, {seconds} seconds"
    datatable = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    datatable.add_column("detail", justify="left", style="green", no_wrap=True)
    datatable.add_column("value", justify="left")
    datatable.add_row("Runtime:", time_text)
    if summary:
        datatable.add_row("Summary: ", f"{outdir}/{summary}")
    console.rule("[bold]Workflow Finished!", style="green")
    console.print(datatable)

def print_onerror(logfile: str) -> None:
    """Print a red panel with error text. To be used in place of onerror: inside a snakefile. Expects the erroring rule printed after it."""
    console.rule("[bold]Workflow Error", style = "red")
    console.print(f"The workflow stopped because of an error. Full workflow log:\n[bold]{logfile}[/]")
    console.rule("[bold]Where Error Occurred", style = "red")

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