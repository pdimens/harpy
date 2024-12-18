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

def print_setup_error(exitcode):
    """Print a red panel with snakefile or conda/singularity error text"""
    if exitcode == 1:
        errortext = "Something is wrong with the Snakefile for this workflow. If you manually edited the Snakefile, see the error below for troubleshooting. If you didn't, it's probably a bug (oops!) and you should submit an issue on GitHub: [bold]https://github.com/pdimens/harpy/issues"
        errortype = "Snakefile Error"
    else:
        errortext = "There was an issue creating the software environment(s) necessary to run this workflow. If you manually edited the conda dependencies in [blue]/workflows/envs[/blue], see the error below for troubleshooting. If you didn't, it might be a bug or related to how your system is setup for Conda or Singularity environments and you should submit an issue on GitHub: [bold]https://github.com/pdimens/harpy/issues"
        errortype = "Software Environment Error"

    rprint(
        Panel(
            errortext,
            title = "[bold]" + errortype,
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
            title = "[dim]Notice",
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

def print_onsuccess(outdir, summary = None):
    """Print a green panel with success text. To be used in place of onsuccess: inside a snakefile"""
    text = f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]"
    if summary:
        text += f" and an outline of the workflow in [bold]{summary}[/bold]"
    rprint(
        Panel(
            text,
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