import sys
from rich import print
from rich.panel import Panel

## define some rich print functions for less redundancy
def print_onstart(text, title):
    """Print a panel of info on workflow run"""
    click.echo("")
    print(Panel(text, title = f"[bold]Harpy {title}", title_align = "left", border_style = "light_steel_blue", subtitle= "Running Workflow"), file = sys.stderr)

def print_error(errortext):
    """Print a yellow panel with error text"""
    print(Panel(errortext, title = "[bold]Error", title_align = "left", border_style = "yellow"), file = sys.stderr)

def print_solution(solutiontext):
    """Print a blue panel with solution text"""
    print(Panel(solutiontext, title = "[bold]Solution", title_align = "left", border_style = "blue"), file = sys.stderr)

def print_solution_with_culprits(solutiontext, culprittext):
    """Print a blue panel with solution text and the list of offenders below it"""
    print(Panel(solutiontext, title = "[bold]Solution", title_align = "left", subtitle = culprittext, border_style = "blue"), file = sys.stderr)

def print_notice(noticetext):
    """Print a white panel with information text text"""
    print(Panel(noticetext, title = "Notice", title_align = "left", border_style = "dim"), file = sys.stderr)