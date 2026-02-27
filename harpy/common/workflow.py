"""Contains the WorkFlow class"""

from datetime import datetime
import importlib.resources as resources
import os
import shutil
import sys
import time as _time
import yaml
from harpy.common.environments import HarpyEnvs
from harpy.common.file_ops import filepath, last_sm_log, purge_empty_logs
from harpy.common.printing import HarpyPrint 
from harpy.common.launch import LaunchSnakemake
from harpy.common.summaries import Summary

class Workflow():
    '''
    The container for workflow parameters. Set inputdir = True to create a workflow/input directory
    '''
    def __init__(self, name, snakefile, outdir, container, clean, quiet, inputdir = False, no_validation: bool = False):
        os.makedirs(
            os.path.join(outdir, 'workflow') if not inputdir else os.path.join(outdir, 'workflow', 'input'),
            exist_ok = True
        )
        self.name: str = name
        self.output_directory: str = outdir
        self.workflow_directory = os.path.join(outdir, 'workflow')
        self.quiet: int = quiet
        self.print = HarpyPrint(quiet)

        self.snakemake_cmd_absolute: str = ""
        self.snakemake_cmd_relative: str = ""
        self.snakefile: str = snakefile

        self.notebook_files: list[str] = []
        self.notebooks: dict = {}
        self.scripts: list[str] = []
        self.inputs: dict = {}
        self.linkedreads: dict = {}
        self.parameters: dict = {}
        self.config: dict = {"Workflow": {}, "Parameters": {}, "Inputs": {}}
        self.profile: dict = {}
        self.hpc: str = ""
        self.clean: str = clean if clean else ""
        self.container: bool = container
        self.conda: list[str] = []

        self.info: dict = {}
        self.start_time: datetime = datetime.now()
        self.summary: str = name.replace("_",".").replace(" ",".") + ".summary"
        self.summary_text: str = ""

        if self.quiet == 0 and not no_validation:
            self.print.console.rule("[bold]Validations", style = "dim magenta")

    def param(self, value, name: str):
        """
        Add the key `name` to self.parameters. if ':' is in the name, then add the left side as
        a key and the right side as a key within that new dict. e.g. `tigmint:mismatch` will add
        `tigmint: {'mismatch: value}` to self.parameters (it will append to the `tigmint` dict if it
        already exists)
        """
        if ":" in name:
            key,subkey = name.split(":")
            if key.strip() not in self.parameters:
                self.parameters[key.strip()] = {}
            self.parameters[key.strip()][subkey.strip()] = value
        else:
            self.parameters[name] = value

    def input(self, value, name: str = "_list"):
        """
        Add the key `name` to self.inputs. if ':' is in the name, then add the left side as
        a key and the right side as a key within that new dict. e.g. `vcf:something` will add
        `vcf: {'something: value}` to self.inputs (it will append to the `tigmint` dict if it
        already exists). If no name is provided, it is assigned to "_list", which will be written
        as a list in the workflow.yaml file.
        """
        if ":" in name:
            key,subkey = name.split(":")
            if key.strip() not in self.inputs:
                self.inputs[key.strip()] = {}
            self.inputs[key.strip()][subkey.strip()] = value
        else:
            self.inputs[name] = value

    def setup_snakemake(self, threads: int, hpc: str|None = None, sm_extra: str|None = None, notemp: bool = False):
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
            formatted_path = os.path.join(*badpath)
            self.print.error(
                "unsupported path name",
                "The path to the output directory includes one or more directories with a space in the name, which is guaranteed to cause errors.",
                f"Rename the path such that there are no spaces in the name:\n{formatted_path}"
            )
        self.profile = {
            "cores" : threads,
            "rerun-incomplete": True,
            "show-failed-logs": True,
            "rerun-triggers": ["mtime", "params"],
            "quiet": "reason",
            "scheduler": "greedy",
            "nolock": True,
            "software-deployment-method": "conda" if not self.container else "apptainer",
            "conda-prefix": filepath("./.environments"),
            "conda-cleanup-pkgs": "cache",
            "apptainer-prefix": filepath("./.environments"),
            "directory": self.output_directory
        }
        _command = ["snakemake", "--snakefile", os.path.join(self.workflow_directory, "workflow.smk")]
        _command += ["--configfile", os.path.join(self.workflow_directory, "workflow.yaml"), "--profile", self.workflow_directory]
        workdir_rel = os.path.relpath(self.workflow_directory)
        _command_rel = ["snakemake", "--snakefile", os.path.join(workdir_rel, "workflow.smk")]
        _command_rel += ["--configfile", os.path.join(workdir_rel, "workflow.yaml"), "--profile", workdir_rel]
        if hpc:
            self.hpc = hpc
            hpc_dir = os.path.join(self.workflow_directory, "hpc")
            self.profile["workflow-profile"] = hpc_dir
        if notemp:
            _command.append("--no-temp")
            _command_rel.append("--no-temp")
        if sm_extra:
            _command.append(sm_extra)
            _command_rel.append(sm_extra)

        self.snakemake_cmd_absolute = " ".join(_command)
        self.snakemake_cmd_relative = " ".join(_command_rel)

    def fetch_notebooks(self) -> None:
        """
        Copy any files in self.notebook_files into workdir/report
        """
        os.makedirs(self.workflow_directory, exist_ok= True)
        
        for target in self.notebook_files:
            dest_file = os.path.join(self.workflow_directory, target)
            source_file = resources.files("harpy.report") / target
            try:
                with resources.as_file(source_file) as _source:
                    shutil.copy2(_source, dest_file)
            except (FileNotFoundError, KeyError):
                self.print.error(
                    "report notebook missing",
                    f"The required report notebook [blue bold]{target}[/] was not found within the Harpy installation.",
                    "There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration."
                )

    def fetch_snakefile(self):
        """
        Copy the workflow snakefile into workflow/workflow.smk
        """
        os.makedirs(self.workflow_directory, exist_ok= True)
        dest_file = os.path.join(self.workflow_directory,"workflow.smk")
        source_file = resources.files("harpy.snakefiles") / self.snakefile
        try:
            with resources.as_file(source_file) as _source:
                shutil.copy2(_source, dest_file)
        except (FileNotFoundError, KeyError):
            self.print.error(
                "snakefile missing",
                f"The required snakefile [blue bold]{self.snakefile}[/] was not found in the Harpy installation.",
                "There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be an issue with your conda/mamba environment or configuration."
            )

    def fetch_scripts(self) -> None:
        """
        Copy any files in self.scripts into workdir/scripts
        """
        for target in self.scripts:
            dest_file = os.path.join(self.workflow_directory, "scripts", target)
            source_file = resources.files("harpy.scripts") / target
            try:
                with resources.as_file(source_file) as _source:
                    shutil.copy2(_source, dest_file)
            except (FileNotFoundError, KeyError):
                self.print.error(
                    "script missing",
                    f"The required script [blue bold]{target}[/] was not found in the Harpy installation.",
                    "There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be an issue with your conda/mamba environment or configuration."
                )

    def fetch_hpc(self):
        """If self.hpc exists, copy it into `workflow/hpc/config.yaml`"""
        if not self.hpc:
            return
        hpc_dest = os.path.join(self.workflow_directory, "hpc")
        os.makedirs(hpc_dest, exist_ok=True)
        shutil.copy2(self.hpc, os.path.join(hpc_dest, "config.yaml"))

    def write_snakemake_profile(self):
        """Writes the Snakemake profile to a file. The profile is expected to be a dict"""
        with open(os.path.join(self.workflow_directory, 'config.yaml'), "w", encoding="utf-8") as sm_config:
            yaml.dump(self.profile, sm_config, sort_keys=False, width=float('inf'))

    def write_workflow_config(self, writefile: bool = True) -> None:
        """
        Formats the workflow configurations into the self.config dict. Writes a workflow.yaml
        file to workdir to use with --configfile if `writefile = True` (default)
        """
        self.config["Workflow"]["name"] = self.name
        if self.linkedreads:
            self.config["Workflow"]["linkedreads"] = self.linkedreads
        if self.notebooks:
            self.config["Workflow"]["reports"] = self.notebooks
        self.config["Workflow"]["snakemake"] = {
            "absolute": self.snakemake_cmd_absolute,
            "relative": self.snakemake_cmd_relative,
            "conda-envs": self.conda
        }
        if self.parameters:
            self.config["Parameters"] = self.parameters
        if "_list" in self.inputs:
            self.config["Inputs"] = self.inputs["_list"]
        else:
            self.config["Inputs"] = self.inputs
        if writefile:
            with open(os.path.join(self.workflow_directory, 'workflow.yaml'), "w", encoding="utf-8") as config:
                yaml.dump(self.config, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    def purge_empty_logs(self):
        purge_empty_logs(self.output_directory)

    def time_elapsed(self) -> str:
        elapsed_time = datetime.now() - self.start_time
        days = elapsed_time.days
        seconds = elapsed_time.seconds
        hours = seconds // 3600
        minutes = (seconds % 3600) // 60
        seconds = seconds % 60
        report_times = []
        for i,j in zip(
            ["day","hour","minute","second"],
            [days,hours,minutes,seconds]
        ):
            if j:
                _i = f"{i}s" if j > 1 else i             
                report_times.append(f"{j} {_i}")
        return ", ".join(report_times)

    def print_onstart(self):
        """Print a panel of info on workflow run. """
        if self.quiet == 2:
            return
        table = self.print.table()
        table.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
        table.add_column("value", justify="left")
        table.add_row("Start:", _time.strftime('%d %b %Y [dim]@[/] %H:%M'))
        for k,v in self.info.items():
                table.add_row(f"{k}:", f"{v}")
        self.print.console.print("")
        self.print.console.rule("[bold]harpy " + self.name.replace("_", " "), style = "light_steel_blue")
        self.print.console.print(table)

    def print_onsuccess(self):
        """Print a green panel with success text. To be used in place of onsuccess: inside a snakefile"""
        if self.quiet == 2:
            return
        _relpath = os.path.relpath(self.output_directory)
        time_text = self.time_elapsed()
        datatable = self.print.table()
        datatable.add_column("detail", justify="left", style="green", no_wrap=True)
        datatable.add_column("value", justify="left")
        datatable.add_row("End:", _time.strftime('%d %b %Y [dim]@[/] %H:%M'))
        datatable.add_row("Duration:", time_text)
        if self.summary:
            datatable.add_row("Summary: ", os.path.join(_relpath, "workflow", os.path.basename(self.summary)))
        _smlog = last_sm_log(self.output_directory)
        if _smlog:
            datatable.add_row("Workflow Log:", os.path.join(_relpath, ".snakemake", 'log', _smlog))
        self.print.rule("[bold]Workflow Finished[/]", style="green")
        self.print.print(datatable)

    def initialize(self, setup: bool = False):
        """Using the configurations, create all necessary folders and files. Launches the workflow if `setup` = False"""
        self.write_workflow_config()
        self.write_snakemake_profile()
        if not self.container:
            HarpyEnvs().write_recipes(self.output_directory, self.conda)
        self.fetch_snakefile()
        self.fetch_scripts()
        self.fetch_notebooks()
        self.fetch_hpc()
        self.print_onstart()
        if not setup:
            self.launch()
        elif self.quiet < 2:
            self.print.rule("[dim bold]Setup Complete", style="dim")

    def launch(self, absolute:bool = False):
        """Launch Snakemake as a monitored subprocess"""
        cmd = self.snakemake_cmd_absolute if absolute else self.snakemake_cmd_relative
        sm = LaunchSnakemake(cmd, self.output_directory, self.quiet, self.print)

        if self.clean:
            self.print.rule("[dim]Cleaning Output", style = "dim")
            for i,j in zip(["w", "s", "l"], ["workflow", ".snakemake", "logs"]):
                if i in self.clean.lower():
                    self.print.log(f"Removing: [blue]{j}/[/]")
                    shutil.rmtree(os.path.join(self.output_directory, j), ignore_errors=True)
        if sm.exitcode == 0:
            with open(os.path.join(self.output_directory, "workflow", f"{self.name.replace('_','.')}.summary"), "w") as f_out: 
                f_out.write(Summary(self.config).get_text())
            self.print_onsuccess()
        else:
            sys.exit(1)
