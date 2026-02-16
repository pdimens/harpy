"""Module of pretty-printing for errors and prompts"""

import time as _time
import os
import re
import sys
from rich.console import Console, RenderableType
from rich import box
from rich.syntax import Syntax
from rich.table import Table
from rich.panel import Panel
from rich.theme import Theme
from harpy.common.version import VERSION

class HarpyPrint():
    def __init__(self, quiet: bool = False):
        '''Instantiate a `HarpyPrint` object configured with a `quiet` setting that applies to printing logs and validations'''
        self.console = Console(
            stderr=True,
            log_path=False,
            theme = Theme({"log.time": "dim magenta"})
        )
        self.quiet = quiet
        self.rule = self.console.rule
        self.print = self.console.print
        self.file = self.console.file
        self.status = self.console.status

    def table(self, title = None, caption = None, tstyle = None, cstyle = None):
        '''
        Insantiate a generic but standardized table style for harpy output
        `tstyle` and `cstyle` refer to the title and caption styles, respectively.
        '''
        return Table(
            title = title,
            caption = caption,
            show_header=False,
            pad_edge=False,
            show_edge=False,
            padding = (0,0),
            box=box.SIMPLE,
            title_style="default" if not tstyle else tstyle,
            caption_style="dim" if not cstyle else cstyle
        )

    def error(
        self,
        errortitle: str,
        errortext: RenderableType|str,
        solutiontext: RenderableType|str|None = None,
        offendertitle: str|None = None,
        offenders: RenderableType|list|str|None = None,
        _exit: bool = True
        ) -> None:
        """
        Print a yellow panel with error text to stderr, exits with error code 1 by default, but exit
        can be disabled with _exit = False. Prints a blue-panel solution if `solutiontext` is provided,
        prints `offenders` after solution if proivded.    
        """
        self.console.print(
            Panel(
                errortext,
                title = f"[bold]Error: {errortitle}",
                title_align = "left",
                border_style = "yellow",
                width = 75
                )
        )
        if solutiontext:
            self.console.print(
                Panel(solutiontext,
                    title = "[bold]Solution",
                    title_align = "left",
                    border_style = "blue",
                    subtitle= offendertitle,
                    width = 75
                )
            )
            if offenders:
                offenders = offenders if isinstance(offenders, str) else "\n".join([str(i) for i in offenders])
                self.console.print(
                    offenders,
                    sep = "\n",
                    highlight= False
                )
        if _exit:
            sys.exit(1)


    def notice(self, noticetext: str|RenderableType) -> None:
        """Print a basic panel with information text to stderr"""
        self.console.print(
            Panel(
                noticetext,
                title = "[dim]Notice",
                title_align = "left",
                border_style = "dim",
                width = 75
                )
        )

    def onstart(self, text: str, title: str) -> None:
        """Print a panel of info on workflow run to stderr"""
        self.console.print("")
        self.console.rule(f"[bold]harpy {title}", style = "light_steel_blue")
        self.console.print(text)

    def setup_error(self, exitcode: int) -> None:
        """Print a red panel with snakefile or conda/singularity error text to stderr"""
        if exitcode == 1:
            errortype = "Snakefile Error"
            errortext = "Something is wrong with the Snakefile for this workflow. If you manually edited the Snakefile, see the error below for troubleshooting. If you didn't, it's probably a bug (oops!) and you should submit an issue on GitHub: [bold]https://github.com/pdimens/harpy/issues"
        else:
            errortype = "Software Environment Error"
            errortext = "There was an issue creating the software environment necessary to run this workflow. If you manually edited the conda dependencies in [blue]/workflows/envs[/], see the error below for troubleshooting. If you didn't, it might be a bug or related to how your system is setup for Conda or Apptainer environments and you should submit an issue on GitHub: [bold]https://github.com/pdimens/harpy/issues"
            # Check if this is the `base` conda environment
            current_env = os.environ.get('CONDA_DEFAULT_ENV')
            if current_env == 'base':
                errortext += "\n[yellow]Notice:[/] Harpy detected that you're in the [blue]base[/] conda environment; conda recommends against installing anything into the [base]base[/] environment."
            # Check the channel priority setting
            try:
                from conda.base.context import context
                if context.channel_priority == "strict":
                    errortext += "\n[yellow]Notice:[/] Your conda channel priority is configured as [yellow]strict[/], which can sometimes cause issues with Snakemake creating conda environments. Ignore this detail if you are using [blue]--container[/]."
            except ModuleNotFoundError:
                pass
        self.console.rule(f"[bold]{errortype}[/][default dim]", style = "red")
        self.console.print("[red]Time:[/] " + _time.strftime('%d %b %Y [dim]@[/] %H:%M'), highlight=False)
        self.console.print(f"[red]Harpy Version:[/] {VERSION}", highlight=False)
        self.console.print(errortext)
        self.console.rule("[bold]Error Reported by Snakemake", style = "red")

    def on_error(self, logfile: str, time) -> None:
        """
        Print a red panel with error text to stderr. To be used in place of onerror: inside a snakefile. Expects the erroring rule printed after it.
        time must be of class datetime
        """
        days = time.days
        seconds = time.seconds
        hours = seconds // 3600
        minutes = (seconds % 3600) // 60
        seconds = seconds % 60
        time_text = f"{days} days, {hours} hours, {minutes} minutes, {seconds} seconds"
        datatable = self.table()
        datatable.add_column("detail", justify="left", style="red", no_wrap=True)
        datatable.add_column("value", justify="left")
        datatable.add_row("Harpy Version:", f"{VERSION}")
        datatable.add_row("Time:", _time.strftime('%d %b %Y @ %H:%M'))
        datatable.add_row("Duration:", time_text)
        datatable.add_row("Workflow Log: ", os.path.relpath(logfile))
        self.console.rule("[bold]Workflow Error[/]", style = "red")
        self.console.print(datatable)
        self.console.print("The workflow stopped due to an error. See the information Snakemake reported below.")
        self.console.rule("[bold]Source of Error", style = "red")

    def shell(self, text) -> None:
        """
        Prints the input text string as syntax-highlighted SHELL code to stderr 
        """
        _table = self.table()
        _table.add_column("Lpadding", justify="left")
        _table.add_column("shell", justify="left")
        _table.add_column("Rpadding", justify="left")

        text = re.sub(r' {2,}|\t+', '  ', text)
        cmd = Syntax(text, lexer = "bash", tab_size=2, word_wrap=True, padding=1, dedent=True, theme = "paraiso-dark")
        _table.add_row("  ", cmd, "  ")
        self.console.print(_table)

    def validation(self, success: bool) -> None:
        '''
        If not `quiet`, print either a red x or green check following a validation log message
        '''
        if not self.quiet:
            if success:
                self.console.print("[green]🗸[/]")
            else:
                self.console.print("[red]𐄂[/]")

    def log(self, text, newline:bool = False):
        '''Print a rich-style log with the time in magenta and text in default'''
        _now =  _time.strftime(r'[dim magenta]\[%H:%M:%S][/]')
        if not self.quiet:
            self.console.print(_now, text, highlight=False, end = "\n" if newline else " ")