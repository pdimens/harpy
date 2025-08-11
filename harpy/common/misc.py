"""Module with helper function to set up Harpy workflows"""

import os
import sys
import glob
import gzip
import shutil
import rich_click as click
from rich.markdown import Markdown
from pathlib import Path
import importlib.resources as resources
from rich.live import Live
from rich.panel import Panel
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, SpinnerColumn, TaskProgressColumn
from .printing import print_error, print_notice, print_solution, CONSOLE

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

def harpy_pulsebar(quiet: int, desc_text: str, stderr: bool = False) -> Progress:
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

def fetch_snakefile(workdir: str, target: str) -> None:
    """
    Retrieve the target harpy rule and write it into the workdir as workflow.smk
    """
    os.makedirs(workdir, exist_ok= True)
    dest_file = os.path.join(workdir,"workflow.smk")
    source_file = resources.files("harpy.snakefiles") / target
    try:
        with resources.as_file(source_file) as _source:
            shutil.copy2(_source, dest_file)
    except (FileNotFoundError, KeyError):
        print_error("snakefile missing", f"The required snakefile [blue bold]{target}[/] was not found in the Harpy installation.", False)
        print_solution("There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration.")

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

def is_conda_package_installed(package_name):
    """
    Return 0 if package exists
    Return 1 if conda is not found but CONDA_PREFIX is defined, implying it's mamba
    Return 2 if package not found in current conda env
    Return 3 if error
    Return 4 if CONDA env not detected
    """
    if "CONDA_PREFIX" in os.environ:
        try:
            from conda.core.prefix_data import PrefixData
            from conda.base.context import context
        except ModuleNotFoundError:
            # CONDA_PREFIX is there but conda itself isnt: likely mamba
            return 1
        try:
            # Get the current environment prefix
            prefix = context.active_prefix or context.root_prefix

            # Get prefix data for current environment
            prefix_data = PrefixData(prefix)

            # Check if package exists in the prefix
            if any(package_name in record.name for record in prefix_data.iter_records()):
                return 0
            else:
                return 2
        except Exception:
            return 3
    else:
        return 4

def is_pip_package_installed(package_name):
    from importlib.metadata import PackageNotFoundError
    import importlib.metadata
    try:
        importlib.metadata.version(package_name)
        return True
    except PackageNotFoundError:
        return False

def is_in_pixi_shell():
    # Check for pixi-specific environment variables
    pixi_indicators = [
        'PIXI_PROJECT_ROOT',
        'PIXI_PROJECT_NAME',
        'PIXI_PROJECT_MANIFEST',
        'PIXI_ENVIRONMENT_NAME'
    ]

    if any(var in os.environ for var in pixi_indicators):
        return True
    return False

def package_absent(pkg: str, executor: bool = True):
    """helper function to search for a package in the active conda environment"""
    if executor:
        out_text = "Using this scheduler requires installing a Snakemake plugin which wasn't detected in this environment. "
    else:
        out_text = "Using the `--container` option requires `apptainer`, which wasn't detected in this environment. "

    out_text += "It can be installed with:"

    # check for conda/mamba
    if shutil.which("conda") or shutil.which("mamba"):
        conda_check = is_conda_package_installed(pkg)
        if conda_check == 0:
            return False
        elif conda_check == 2:
            if is_in_pixi_shell():
                out_text += f"\n\n```bash\npixi add {pkg}\n```"
            else:
                out_text += f"\n\n```bash\nconda install bioconda::{pkg}\n```"
        elif conda_check == 1:
            out_text += f"\n\n```bash\nmamba install bioconda::{pkg}\n```"
        if conda_check in [1,2]:
            if executor:
                print_notice(Markdown(out_text))
            else:
                print_error("missing required package", Markdown(out_text))
            return True

    if not is_pip_package_installed(pkg):
        out_text += f"\n\n```bash\npip install {pkg}\n```"
        if executor:
            print_notice(Markdown(out_text))
        else:
            print_error("missing required package", Markdown(out_text))
        return True
    
    return False

def container_ok(ctx, param, value):
    if value:
        if os.sys.platform != 'linux':
            raise click.BadParameter(
                "Snakemake uses Singularity to manage containers, which is only available for Linux systems.", ctx, param
            )
        if shutil.which("apptainer"):
            return value
        else:
            raise click.BadParameter(
                "Container software management requires apptainer, which wasn't detected in this environment.", ctx, param
            )
