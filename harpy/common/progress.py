"""
Functions related to Harpy's progressbars
"""

from rich.live import Live
from rich.panel import Panel
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, SpinnerColumn, TaskProgressColumn
from rich.text import Text
from harpy.common.printing import CONSOLE
import time

class PausableTimeElapsedColumn(TimeElapsedColumn):
    """Custom time elapsed column that supports pausing and resuming."""
    
    def __init__(self):
        super().__init__()
        self.pause_adjustments = {}  # task_id -> total paused time
        self.pause_start_times = {}  # task_id -> when pause started
    
    def pause(self, task_id):
        """Start pausing the timer for a task."""
        self.pause_start_times[task_id] = time.monotonic()
    
    def resume(self, task_id):
        """Resume the timer for a task."""
        if task_id in self.pause_start_times:
            pause_duration = time.monotonic() - self.pause_start_times[task_id]
            self.pause_adjustments[task_id] = self.pause_adjustments.get(task_id, 0) + pause_duration
            del self.pause_start_times[task_id]
    
    def render(self, task):
        """Render the elapsed time, accounting for pauses."""
        elapsed = task.elapsed
        
        # subtract any paused time
        if task.id in self.pause_adjustments:
            elapsed -= self.pause_adjustments[task.id]

        # if currently paused, also subtract time since pause started
        if task.id in self.pause_start_times:
            elapsed -= (time.monotonic() - self.pause_start_times[task.id])

        # don't go negative
        elapsed = max(0, elapsed)
        
        # Format the time
        minutes, seconds = divmod(int(elapsed), 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)
        
        if days:
            return Text(f"{days:d} days, {hours:d} hours", style="progress.elapsed")
        else:
            return Text(f"{hours:d}:{minutes:02d}:{seconds:02d}", style="progress.elapsed")

def harpy_progresspanel(progressbar: Progress, title: str|None = None, quiet: int = 0):
    """Returns a nicely formatted live-panel with the progress bar in it"""
            #progressbar if quiet != 2 else None,
    return Live(
        Panel(
            progressbar, title = title, border_style="dim"
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
        TextColumn("{task.fields[active]}", style="default"),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=None, complete_style="yellow", finished_style="dim blue"),
        TaskProgressColumn("[progress.remaining]{task.completed}/{task.total}") if quiet == 0 else TaskProgressColumn(),
        PausableTimeElapsedColumn(),
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