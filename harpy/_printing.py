"""Module of pretty-printing for errors and prompts"""

import sys
from rich import print as rprint
from rich.panel import Panel

def print_error(errortitle, errortext):
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

def print_snakefile_error(errortext):
    """Print a red panel with snakefile error text"""
    rprint(
        Panel(
            errortext,
            title = "[bold]Snakefile Error",
            title_align = "left",
            subtitle = "The error reported by Snakemake",
            border_style = "red",
            width = 75
            ),
        file = sys.stderr
    )

def print_solution(solutiontext):
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

def print_solution_with_culprits(solutiontext, culprittext):
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

def print_notice(noticetext):
    """Print a white panel with information text text"""
    rprint(
        Panel(
            noticetext,
            title = "Notice",
            title_align = "left",
            border_style = "dim",
            width = 75
            ),
        file = sys.stderr
    )

def print_onstart(text, title):
    """Print a panel of info on workflow run"""
    rprint("")
    rprint(
        Panel(
            text,
            title = f"[bold]harpy {title}",
            title_align = "center",
            border_style = "light_steel_blue",
            width = 75
            ),
        file = sys.stderr
    )

def print_onsuccess(outdir):
    """Print a green panel with success text. To be used in place of onsuccess: inside a snakefile"""
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]",
            subtitle = "[bold]success!",
            title_align = "left",
            border_style = "green",
            width=75
            ),
        file = sys.stderr
    )

def print_onerror(logfile):
    """Print a red panel with error text. To be used in place of onsuccess: inside a snakefile. Expects the erroring rule printed after it."""
    rprint(
        Panel(
            f"The workflow terminated from an error. See the full log for more info:\n[bold]{logfile}[/bold]",
            title = "[bold]workflow error",
            title_align = "left",
            border_style = "red",
            subtitle = "the step causing this error:",
            width=75,
            ),
        file = sys.stderr
    )