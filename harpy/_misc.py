"""Module with helper function to set up Harpy workflows"""

import os
import sys
import glob
import gzip
import shutil
import yaml
import urllib.request
from datetime import datetime
from pathlib import Path
from importlib_resources import files
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, SpinnerColumn, TaskProgressColumn, MofNCompleteColumn
import harpy.scripts
import harpy.reports
import harpy.snakefiles
from ._printing import print_error, print_solution

def harpy_progressbar(quiet: int) -> Progress:
    """
    The pre-configured transient progress bar that workflows and validations use
    """
    return Progress(
        SpinnerColumn(spinner_name = "dots12", style = "blue dim"),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(complete_style="yellow", finished_style="blue"),
        TaskProgressColumn("[progress.remaining]{task.completed}/{task.total}") if quiet == 0 else TaskProgressColumn(),
        TimeElapsedColumn(),
        transient = True,
        auto_refresh = True,
        disable = quiet == 2
    )

def harpy_pulsebar(quiet: int, desc_text: str) -> Progress:
    """
    The pre-configured transient pulsing progress bar that workflows use, typically for
    installing the software dependencies/container
    """
    return Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width= 70 - len(desc_text), pulse_style = "grey46"),
        TimeElapsedColumn(),
        auto_refresh = True,
        transient = True,
        disable = quiet == 2
    )

def symlink(original: str, destination: str) -> None:
    """Create a symbolic link from original -> destination if the destination doesn't already exist."""
    if not (Path(destination).is_symlink() or Path(destination).exists()):
        Path(destination).symlink_to(Path(original).resolve())

def fetch_script(workdir: str, target: str) -> None:
    """
    Retrieve the target harpy script and write it into workdir/scripts
    """
    os.makedirs(f"{workdir}/scripts/", exist_ok= True)
    with open(f"{workdir}/scripts/{target}", "w", encoding="utf-8") as f:
        if os.path.isfile(files(harpy.scripts).joinpath(target)):
            f.write(files(harpy.scripts).joinpath(target).read_text())
        else:
            print_error("script missing", f"Bundled script [blue bold]{target}[/] was not found in the Harpy installation.")
            print_solution("There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration.")
            sys.exit(1)

def fetch_rule(workdir: str, target: str) -> None:
    """
    Retrieve the target harpy rule and write it into the workdir
    """
    os.makedirs(f"{workdir}/", exist_ok= True)
    with open(f"{workdir}/{target}", "w", encoding="utf-8") as f:
        if os.path.isfile(files(harpy.snakefiles).joinpath(target)):
            f.write(files(harpy.snakefiles).joinpath(target).read_text())
        else:
            print_error("snakefile missing", f"The required snakefile [blue bold]{target}[/] was not found in the Harpy installation.")
            print_solution("There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration.")
            sys.exit(1)

def fetch_report(workdir: str, target: str) -> None:
    """
    Retrieve the target harpy report and write it into workdir/report
    """
    os.makedirs(f"{workdir}/report/", exist_ok= True)
    with open(f"{workdir}/report/{target}", "w", encoding="utf-8") as f:
        if os.path.isfile(files(harpy.reports).joinpath(target)):
            f.write(files(harpy.reports).joinpath(target).read_text())
        else:
            print_error("report script missing", f"The required report script [blue bold]{target}[/] was not found within the Harpy installation.")
            print_solution("There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration.")
            sys.exit(1)
    # pull yaml config file from github, use local if fails
    try:
        yaml = "https://github.com/pdimens/harpy/raw/refs/heads/main/harpy/reports/_quarto.yml"
        destination = f"{workdir}/report/_quarto.yml"
        with urllib.request.urlopen(yaml) as response, open(destination, 'w') as yaml_out:
            yaml_out.write(response.read().decode("utf-8"))
    except:
        with open(f"{workdir}/report/_quarto.yml", "w", encoding="utf-8") as yml:
            if os.path.isfile(files(harpy.reports).joinpath("_quarto.yml")):
                yml.write(files(harpy.reports).joinpath("_quarto.yml").read_text())
            else:
                print_error("report configuration missing", f"The required quarto configuration could not be downloaded from the Harpy repository, nor found in the local file [blue bold]_quarto.yml[/] that comes with a Harpy installation.")
                print_solution("There may be an issue with your internet connection or Harpy installation, that latter of which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration.")
                sys.exit(1)
    # same for the scss file
    try:
        scss = "https://github.com/pdimens/harpy/raw/refs/heads/main/harpy/reports/_harpy.scss"
        destination = f"{workdir}/report/_harpy.scss"
        with urllib.request.urlopen(scss) as response, open(destination, 'w') as scss_out:
            scss_out.write(response.read().decode("utf-8"))
    except:
        with open(f"{workdir}/report/_harpy.scss", "w", encoding="utf-8") as yml:
            if os.path.isfile(files(harpy.reports).joinpath("_harpy.scss")):
                yml.write(files(harpy.reports).joinpath("_harpy.scss").read_text())
            else:
                print_error("report configuration missing", f"The required quarto configuration could not be downloaded from the Harpy repository, nor found in the local file [blue bold]_harpy.scss[/] that comes with a Harpy installation.")
                print_solution("There may be an issue with your internet connection or Harpy installation, that latter of which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration.")
                sys.exit(1)

def snakemake_log(outdir: str, workflow: str) -> str:
    """Return a snakemake logfile name. Iterates logfile run number if one exists."""
    attempts = glob.glob(f"{outdir}/logs/snakemake/*.log*")
    if not attempts:
        return f"logs/snakemake/{workflow}.1." + datetime.now().strftime("%d_%m_%Y") + ".log"
    increment = sorted([int(i.split(".")[1]) for i in attempts])[-1] + 1
    return f"logs/snakemake/{workflow}.{increment}." + datetime.now().strftime("%d_%m_%Y") + ".log"

def gzip_file(infile: str) -> None:
    """gzip a file and delete the original, using only python"""
    if os.path.exists(infile):
        with open(infile, 'rb') as f_in, gzip.open(infile + '.gz', 'wb', 6) as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(infile)

def safe_read(file_path: str):
    """returns the proper file opener for reading if a file_path is gzipped"""
    try:
        with gzip.open(file_path, 'rt') as f:
            f.read(10)
        return gzip.open(file_path, 'rt')
    except gzip.BadGzipFile:
        return open(file_path, 'r')

def write_snakemake_config(sdm, outdir):
    """
    Writes a config.yaml file to outdir/workflow to use with --profile. Creates outdir/workflow if it doesnt exist.
    sdm is the software deployment method
    """
    profile = {
        "rerun-incomplete": True,
        "show-failed-logs": True,
        "rerun-triggers": ["mtime", "params"],
        "nolock": True,
        "software-deployment-method": sdm,
        "conda-prefix": Path("./.environments").resolve().as_posix(),
        "conda-cleanup-pkgs": "cache",
        "apptainer-prefix": Path("./.environments").resolve().as_posix(),
        "directory": outdir
    }
    workdir = os.path.join(outdir, "workflow")
    if not os.path.exists(workdir):
        os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, 'config.yaml'), "w", encoding="utf-8") as sm_config:
        yaml.dump(profile, sm_config, sort_keys=False, width=float('inf'))

def write_workflow_config(configs, outdir):
    """
    Writes a workflow.yaml file to workdir to use with --configfile. Creates outdir/workflow if it doesnt exist. Configs
    are expected to be a dict
    """
    workdir = f"{outdir}/workflow"
    if not os.path.exists(workdir):
        os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, f'config.harpy.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))