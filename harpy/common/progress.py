"""progress.py"""
from rich.table import Column
from rich.text import Text
import time
from rich.progress import TimeElapsedColumn
from rich.progress import BarColumn
from rich.progress import TaskProgressColumn
from rich.progress import TextColumn
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress

#TODO MAKE PULSEBAR TRANSIENT
class PanelProgress(Progress):
    def __init__(self, console: Console, quiet: int, title: str | None = None, border_style: str = "dim", transient: bool = False):
        self.quiet = quiet
        self.panel_title = title
        self._border_style = border_style
        super().__init__(
            console=console,
            transient=transient,
            disable=quiet == 2,
            auto_refresh=True,
            refresh_per_second=2,
            expand=True,
        )

    def get_renderables(self):
        yield Panel(
            self.make_tasks_table(self.tasks),
            title=self.panel_title,
            border_style=self._border_style,
            expand=True,
        )

    def start(self):
        if self.quiet != 2:
            super().start()

    def stop(self):
        if self.quiet != 2:
            super().stop()

    def add_task(self, *args, **kwargs):
        return -1 if self.quiet == 2 else super().add_task(*args, **kwargs)

    def update(self, task_id, **kwargs):
        if self.quiet != 2:
            super().update(task_id, **kwargs)

    def bar(self) -> "PanelProgress":
        '''A progressbar for tracking snakemake jobs'''
        self.columns = (
            TextColumn("{task.fields[active]}", style="yellow"),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(bar_width=None, complete_style="yellow", finished_style="dim blue"),
            TaskProgressColumn("{task.completed}/{task.total}", style="blue") if self.quiet == 0 else TaskProgressColumn(style="blue"),
            PausableTimeElapsedColumn(),
        )
        return self

    def basic(self, width = None) -> "PanelProgress":
        '''A simpler progress bar for tracking a singular task'''
        self.columns = (
            TextColumn("[progress.description]{task.description}", table_column=Column(width=width)),
            BarColumn(bar_width=None, complete_style="yellow", finished_style="dim blue"),
            TaskProgressColumn("{task.completed}/{task.total}", style="blue"),
            TimeElapsedColumn(),
        )
        return self


    def pulse(self) -> "PanelProgress":
        '''A pulsing progress bar used for conda/apptainer progress'''
        self.columns = (
            TextColumn("[progress.description]{task.description}"),
            BarColumn(bar_width=None, pulse_style="grey46"),
            TimeElapsedColumn(),
        )
        return self


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
        _style = "yellow"

        # subtract any paused time
        if task.id in self.pause_adjustments:
            elapsed -= self.pause_adjustments[task.id]

        # if currently paused, also subtract time since pause started
        if task.id in self.pause_start_times:
            elapsed -= (time.monotonic() - self.pause_start_times[task.id])
            _style = "dim yellow"

        # don't go negative
        elapsed = max(0, elapsed)

        # Format the time
        minutes, seconds = divmod(int(elapsed), 60)
        hours, minutes = divmod(minutes, 60)
        days, hours = divmod(hours, 24)

        if days:
            _days = "day" if days == 1 else "days"
            _hours = "hour" if hours == 1 else "hours"
            return Text(f"{days:d} {_days}, {hours:d} {_hours}", style = _style)
        else:
            return Text(f"{hours:d}:{minutes:02d}:{seconds:02d}", style = _style)
