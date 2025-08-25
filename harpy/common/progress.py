"""
Functions related to Harpy's progressbars
"""

from rich.live import Live
from rich.panel import Panel
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, SpinnerColumn, TaskProgressColumn
from harpy.common.printing import CONSOLE

def harpy_progresspanel(progressbar: Progress, title: str|None = None, quiet: int = 0):
    """Returns a nicely formatted live-panel with the progress bar in it"""
    return Live(
        Panel(
            progressbar if quiet != 2 else None,
            title = title,
            border_style="dim"
        ) if quiet != 2 else None,
        refresh_per_second=8,
        transient=True,
        console=CONSOLE
    )

def harpy_progressbar(quiet: int) -> Progress:
    """
    The pre-configured transient progress bar that workflows and validations use
    """
    return Progress(
        SpinnerColumn(spinner_name = "dots12", style = "blue dim", finished_text="[dim green]✓"),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=None, complete_style="yellow", finished_style="dim blue"),
        TaskProgressColumn("[progress.remaining]{task.completed}/{task.total}") if quiet == 0 else TaskProgressColumn(),
        TimeElapsedColumn(),
        transient = True,
        auto_refresh = True,
        disable = quiet == 2,
        console= CONSOLE,
        expand=True
    )

def harpy_pulsebar(quiet: int, stderr: bool = False) -> Progress:
    """
    The pre-configured transient pulsing progress bar that workflows use, typically for
    installing the software dependencies/container
    """
    return Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width= None, pulse_style = "grey46"),
        TimeElapsedColumn(),
        auto_refresh = True,
        transient = True,
        disable = quiet == 2,
        console = CONSOLE if stderr else None,
        expand=True
    )