"""Contains the WorkFlow class"""

from datetime import datetime
import glob
import importlib.resources as resources
import os
import shutil
import time as _time
import urllib.request
import yaml
from rich import print as rprint
from rich import box
from rich.table import Table
from harpy.common.conda import create_conda_recipes
from harpy.common.file_ops import filepath, gzip_file, fetch_snakefile, purge_empty_logs
from harpy.common.printing import CONSOLE, print_error
from harpy.common.launch import launch_snakemake

class Workflow():
    '''
    The container for workflow parameters. Set inputdir = True to create a workflow/input directory
    '''
    def __init__(self, name, snakefile, outdir, quiet, inputdir = False):
        creatdir = os.path.join(outdir, 'workflow') if not inputdir else os.path.join(outdir, 'workflow', 'input')
        os.makedirs(creatdir, exist_ok = True)

        self.name = name
        self.output_directory = outdir
        self.workflow_directory = os.path.join(outdir, 'workflow')
        self.snakemake_log = self.snakemake_log(outdir, name)
        self.snakemake_cmd_absolute = None
        self.snakemake_cmd_relative = None
        self.snakefile = snakefile
        self.reports =[]
        self.scripts = []
        self.config = None
        self.profile = None
        self.hpc = None
        self.conda = None
        self.start_text = None
        self.quiet = quiet
        self.start_time = datetime.now()
        self.summary = name.replace("_",".").replace(" ",".") + ".summary"

    def snakemake_log(self, outdir: str, workflow: str) -> str:
        """Return a snakemake logfile name. Iterates logfile run number if one exists."""
        attempts = glob.glob(os.path.join(outdir, "logs", "snakemake", "*.log*"))
        timestamp = datetime.now().strftime("%d_%m_%Y") + ".log"
        if not attempts:
            os.makedirs(os.path.join(outdir, "logs", "snakemake"), exist_ok=True)
            return os.path.join("logs", "snakemake", f"{workflow}.1.{timestamp}")
        increment = sorted([int(i.split(".")[1]) for i in attempts])[-1] + 1
        return os.path.join("logs", "snakemake", f"{workflow}.{increment}.{timestamp}")

    #def setup_snakemake(self, sdm: str, threads: int, hpc: str|None = None, sm_extra: str|None = None):
    def setup_snakemake(self, container: bool, threads: int, hpc: str|None = None, sm_extra: str|None = None):
        """
        Sets up the snakemake command based on hpc, threads, and extra snakemake params.
        """
        badpath = []
        patherr = False
        for i in self.output_directory.split("/"):
            if " " in i:
                patherr = True
                badpath.append(f"[red]{i}[/]")
            else:
                badpath.append(i)
        if patherr:
            formatted_path = os.path.join(badpath)
            print_error(
                "unsupported path name",
                "The path to the output directory includes one or more directories with a space in the name, which is guaranteed to cause errors.",
                f"Rename the path such that there are no spaces in the name:\n{formatted_path}"
            )
        self.profile = {
            "rerun-incomplete": True,
            "show-failed-logs": True,
            "rerun-triggers": ["mtime", "params"],
            "nolock": True,
            "software-deployment-method": "conda" if not container else "conda apptainer",
            "conda-prefix": filepath("./.environments"),
            "conda-cleanup-pkgs": "cache",
            "apptainer-prefix": filepath("./.environments"),
            "directory": self.output_directory
        }
        _command = []
        _command += ["snakemake", "--cores", f"{threads}", "--snakefile", os.path.join(self.workflow_directory, "workflow.smk")]
        _command += ["--configfile", os.path.join(self.workflow_directory, "workflow.yaml"), "--profile", self.workflow_directory]
        workdir_rel = os.path.relpath(self.workflow_directory)
        _command_rel = []
        _command_rel += ["snakemake", "--cores", f"{threads}", "--snakefile", os.path.join(workdir_rel, "workflow.smk")]
        _command_rel += ["--configfile", os.path.join(workdir_rel, "workflow.yaml"), "--profile", workdir_rel]
        if hpc:
            self.hpc = hpc
            hpc_dir = os.path.join(self.workflow_directory, "hpc")
            os.makedirs(hpc_dir, exist_ok=True)
            _command += ["--workflow-profile", hpc_dir]
            _command_rel += ["--workflow-profile", os.path.join(workdir_rel, "hpc")]
        if sm_extra:
            _command.append(sm_extra)
            _command_rel.append(sm_extra)
    
        self.snakemake_cmd_absolute = " ".join(_command)
        self.snakemake_cmd_relative = " ".join(_command_rel)

    def fetch_report(self, target: str) -> None:
        """
        Retrieve the target harpy report and write it into workdir/report
        """
        dest_dir = os.path.join(self.workflow_directory, "report")
        dest_file = os.path.join(dest_dir, target)
        os.makedirs(dest_dir, exist_ok= True)
        source_file = resources.files("harpy.reports") / target
        try:
            with resources.as_file(source_file) as _source:
                shutil.copy2(_source, dest_file)
        except (FileNotFoundError, KeyError):
            print_error(
                "report script missing",
                f"The required report script [blue bold]{target}[/] was not found within the Harpy installation.",
                "There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration."
            )

    def fetch_report_configs(self):
        """
        pull yaml config file from github, use local if fails
        pull yaml config file from GitHub, use local if download fails
        """
        dest_dir = os.path.join(self.workflow_directory, "report")
        destination = os.path.join(dest_dir, "_quarto.yml")
        try:
            _yaml = "https://github.com/pdimens/harpy/raw/refs/heads/main/harpy/reports/_quarto.yml"
            with urllib.request.urlopen(_yaml) as response, open(destination, 'w') as yaml_out:
                yaml_out.write(response.read().decode("utf-8"))
        except (urllib.error.URLError, OSError, IOError):
            source_file = resources.files("harpy.reports") / "_quarto.yml"
            try:
                with resources.as_file(source_file) as _source:
                    shutil.copy2(_source, destination)
            except (FileNotFoundError, KeyError):
                print_error(
                    "report configuration missing",
                    "The required quarto configuration could not be downloaded from the Harpy repository, nor found in the local file [blue bold]_quarto.yml[/] that comes with a Harpy installation.",
                    "There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be an issue with your conda/mamba environment or configuration."
                )

        # same for the scss file
        destination = os.path.join(dest_dir, "_harpy.scss")
        try:
            scss = "https://github.com/pdimens/harpy/raw/refs/heads/main/harpy/reports/_harpy.scss"
            with urllib.request.urlopen(scss) as response, open(destination, 'w') as scss_out:
                scss_out.write(response.read().decode("utf-8"))
        except (urllib.error.URLError, OSError, IOError):
            source_file = resources.files("harpy.reports") / "_harpy.scss"
            try:
                with resources.as_file(source_file) as _source:
                    shutil.copy2(_source, destination)
            except (FileNotFoundError, KeyError):
                print_error(
                    "report configuration missing",
                    "The required quarto configuration could not be downloaded from the Harpy repository, nor found in the local file [blue bold]_harpy.scss[/] that comes with a Harpy installation.",
                    "There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration."
                    )

    def fetch_snakefile(self):
         """
         Retrieve the target harpy rule and write it into the workdir as workflow.smk
         """
         fetch_snakefile(self.workflow_directory, self.snakefile)

    def fetch_script(self, target: str) -> None:
        """
        Retrieve the target harpy script and write it into workdir/scripts
        """
        dest_file = os.path.join(self.workflow_directory, "scripts", target)
        source_file = resources.files("harpy.scripts") / target
        try:
            with resources.as_file(source_file) as _source:
                shutil.copy2(_source, dest_file)
        except (FileNotFoundError, KeyError):
            print_error(
                "snakefile missing",
                f"The required script [blue bold]{target}[/] was not found in the Harpy installation.",
                "There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be an issue with your conda/mamba environment or configuration."
            )

    def fetch_hpc(self):
        hpc_dest = os.path.join(self.workflow_directory, "hpc")
        os.makedirs(hpc_dest, exist_ok=True)
        shutil.copy2(self.hpc, os.path.join(hpc_dest, "config.yaml"))

    def write_snakemake_profile(self):
        """Writes the Snakemake profile to a file. The profile is expected to be a dict"""
        with open(os.path.join(self.workflow_directory, 'config.yaml'), "w", encoding="utf-8") as sm_config:
            yaml.dump(self.profile, sm_config, sort_keys=False, width=float('inf'))

    def write_workflow_config(self) -> None:
        """
        Writes a workflow.yaml file to workdir to use with --configfile. Configs
        are expected to be a dict
        """
        with open(os.path.join(self.workflow_directory, 'workflow.yaml'), "w", encoding="utf-8") as config:
            yaml.dump(self.config, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    def purge_empty_logs(self):
        purge_empty_logs(self.output_directory)

    def initialize(self, setup_only: bool = False):
        """Using the configurations, create all necessary folders and files"""
        self.write_workflow_config()
        self.write_snakemake_profile()
        create_conda_recipes(self.output_directory, self.conda)
        self.fetch_snakefile()
        for i in self.reports:
            self.fetch_report(i)
        for i in self.scripts:
            self.fetch_script(i)
        if self.reports:
            self.fetch_report_configs()
        if self.hpc:
            self.fetch_hpc()
        self.print_onstart()
        if not setup_only:
            self.launch()

    def print_onstart(self):
        """Print a panel of info on workflow run"""
        if self.quiet == 2:
            return
        rprint("")
        CONSOLE.rule("[bold]harpy " + self.name.replace("_", " "), style = "light_steel_blue")
        CONSOLE.print(self.start_text)

    def print_onsuccess(self):
        """Print a green panel with success text. To be used in place of onsuccess: inside a snakefile"""
        if self.quiet == 2:
            return
        elapsed_time = datetime.now() - self.start_time 
        days = elapsed_time.days
        seconds = elapsed_time.seconds
        hours = seconds // 3600
        minutes = (seconds % 3600) // 60
        seconds = seconds % 60
        time_text = f"{days} days, {hours} hours, {minutes} minutes, {seconds} seconds"
        datatable = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
        datatable.add_column("detail", justify="left", style="green", no_wrap=True)
        datatable.add_column("value", justify="left")
        datatable.add_row("Duration:", time_text)
        if self.summary:
            datatable.add_row("Summary: ", f"{os.path.basename(self.output_directory)}/workflow/{self.summary}")
        datatable.add_row("Workflow Log:", f"{os.path.basename(self.output_directory)}/{self.snakemake_log}.gz")
        CONSOLE.rule("[bold]Workflow Finished[/] [dim]" + _time.strftime('%d %b %Y @ %H:%M'), style="green")
        CONSOLE.print(datatable)

    def launch(self, absolute:bool = False):
        """Launch Snakemake as a monitored subprocess"""
        if absolute:
            launch_snakemake(self.snakemake_cmd_absolute, self.workflow_directory, self.output_directory, self.snakemake_log, self.quiet)
        else:
            launch_snakemake(self.snakemake_cmd_relative, self.workflow_directory, self.output_directory, self.snakemake_log, self.quiet)
        self.purge_empty_logs()
        gzip_file(os.path.join(self.output_directory, self.snakemake_log))
        self.print_onsuccess()