"""Module with helper function to set up Harpy workflows"""

import os
import sys
import glob
import gzip
import shutil
from pathlib import Path
import importlib.resources as resources
from rich.live import Live
from rich.panel import Panel
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, SpinnerColumn, TaskProgressColumn
from rich.console import Console
from .printing import print_error, print_solution

_STDERR_CONSOLE = Console(file=sys.stderr)

def harpy_progresspanel(progressbar: Progress, title: str|None = None, quiet: int = 0):
    """Returns a nicely formatted live-panel with the progress bar in it"""
    return Live(
        Panel(
            progressbar if quiet != 2 else None,
            title = title,
            border_style="dim"
        ) if quiet != 2 else None,
        refresh_per_second=8,
        transient=True
    )

def harpy_progressbar(quiet: int) -> Progress:
    """
    The pre-configured transient progress bar that workflows and validations use
    """
    return Progress(
        SpinnerColumn(spinner_name = "dots12", style = "blue dim", finished_text="[dim green]✓"),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(complete_style="yellow", finished_style="dim blue"),
        TaskProgressColumn("[progress.remaining]{task.completed}/{task.total}") if quiet == 0 else TaskProgressColumn(),
        TimeElapsedColumn(),
        transient = True,
        auto_refresh = True,
        disable = quiet == 2
    )

def harpy_pulsebar(quiet: int, desc_text: str, stderr: bool = False) -> Progress:
    """
    The pre-configured transient pulsing progress bar that workflows use, typically for
    installing the software dependencies/container
    """
    return Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width= max(10, 70 - len(desc_text)), pulse_style = "grey46"),
        TimeElapsedColumn(),
        auto_refresh = True,
        transient = True,
        disable = quiet == 2,
        console = _STDERR_CONSOLE if stderr else None
    )

def fetch_snakefile(workdir: str, target: str) -> None:
    """
    Retrieve the target harpy rule and write it into the workdir as workflow.smk
    """
    dest_file = os.path.join(workdir,"workflow.smk")
    os.makedirs(workdir, exist_ok= True)
    source_file = resources.files("harpy.snakefiles") / target
    try:
        with resources.as_file(source_file) as _source:
            shutil.copy2(_source, dest_file)
    except (FileNotFoundError, KeyError):
        print_error("snakefile missing", f"The required snakefile [blue bold]{target}[/] was not found in the Harpy installation.")
        print_solution("There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration.")
        sys.exit(1)

def filepath(infile: str) -> str:
    """returns a posix-formatted absolute path of infile"""
    return Path(infile).resolve().as_posix()

def symlink(original: str, destination: str) -> None:
    """Create a symbolic link from original -> destination if the destination doesn't already exist."""
    if not (Path(destination).is_symlink() or Path(destination).exists()):
        Path(destination).symlink_to(Path(original).resolve())

def gzip_file(infile: str) -> None:
    """gzip a file and delete the original, using only python"""
    if os.path.exists(infile):
        with open(infile, 'rb') as f_in, gzip.open(infile + '.gz', 'wb', 6) as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(infile)

def purge_empty_logs(output_directory):
    """scan target_dir and remove empty files, then scan it again and remove empty directories"""
    for logfile in glob.glob(f"{output_directory}/logs/**/*", recursive = True):
        if os.path.isfile(logfile) and os.path.getsize(logfile) == 0:
            os.remove(logfile)
    for logfile in glob.glob(f"{output_directory}/logs/**/*", recursive = True):
        if os.path.isdir(logfile) and not os.listdir(logfile):
            os.rmdir(logfile)

def safe_read(file_path: str):
    """returns the proper file opener for reading if a file_path is gzipped"""
    try:
        with gzip.open(file_path, 'rt') as f:
            f.read(10)
        return gzip.open(file_path, 'rt')
    except gzip.BadGzipFile:
        return open(file_path, 'r')
