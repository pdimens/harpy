"""Module of pretty-printing for errors and prompts"""

import os
import re
import sys
import time
from contextlib import nullcontext

from rich import box
from rich.console import Console, RenderableType
from rich.live import Live
from rich.markup import escape
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text
from rich.theme import Theme

from harpy import __version__
from .ruleparsing import _Pushback, SnakeRule, format_missing_output_block, _is_separator, format_shell_cmd
from .ruleparsing import _LOGFILE_NOTFOUND_RE, _HEADER_RE, _KEY_RE, _LOGFILE_HEADER_RE

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

class HarpyPrint():
    def __init__(self, quiet: int = 0):
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

    def time_now(self) -> str:
        return time.strftime('%H:%M [dim]on[/] %d %b %Y')

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
        self.print(
            Panel(
                errortext,
                title = f"[bold]Error: {errortitle}",
                title_align = "left",
                border_style = "yellow",
                width = 75
                )
        )
        if solutiontext:
            self.print(
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
                self.print(
                    offenders,
                    sep = "\n",
                    highlight= False
                )
        if _exit:
            sys.exit(1)

    def notice(self, noticetext: str|RenderableType) -> None:
        """Print a basic panel with information text to stderr"""
        self.print(
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
        self.print("")
        self.rule(f"[bold]{title}", style = "light_steel_blue")
        self.print(text)

    def setup_error(self, exitcode: int) -> None:
        """Print a red panel with snakefile or conda/singularity error text to stderr"""
        if exitcode == 1:
            errortype = "Snakefile Error"
            errortext = "Something is wrong with this Snakefile. If you edited it, see the error below for troubleshooting, and if not, please submit an issue on GitHub: [bold]https://github.com/pdimens/harpy/issues"
        else:
            errortype = "Software Environment Error"
            errortext = "There was an issue creating the environment for this workflow. If you manually edited the Conda dependencies in [blue]/workflows/envs[/], see the error below for troubleshooting. Otherwise, it may be a bug or related to your Conda/Apptainer setup—please report it on GitHub: [bold]https://github.com/pdimens/harpy/issues[/]"
            # Check if this is the `base` conda environment
            current_env = os.environ.get('CONDA_DEFAULT_ENV')
            if current_env == 'base':
                errortext += "\n[yellow]Notice:[/] Harpy detected that you're in the [blue]base[/] conda environment; conda recommends against installing anything into the [blue]base[/] environment."
            # Check the channel priority setting
            try:
                from conda.base.context import context
                if context.channel_priority == "strict":
                    errortext += "\n[yellow]Notice:[/] Your conda channel priority is configured as [yellow]strict[/], which can sometimes cause issues with Snakemake creating conda environments. Ignore this detail if you are using [blue]--container[/]."
            except ModuleNotFoundError:
                pass
        self.print("")
        self.rule(f"[bold]{errortype}[/][default dim]", style = "red")
        self.print(f"[red]Harpy:[/] v{__version__}", highlight=False)
        self.print("[red]Time:[/] " + self.time_now(), highlight=False)
        self.print(errortext + "\n")
        self.print("[bold dim]──── ⚠ Error Reported by Snakemake")

    def on_error(self, logfile: str, _time) -> None:
        """
        Print a red panel with error text to stderr. To be used in place of onerror: inside a snakefile. Expects the erroring rule printed after it.
        time must be of class datetime
        """
        #TODO NO LONGER CATCHES SNAKEFILE ERRORS. ADDRESS THAT
        days = _time.days
        seconds = _time.seconds
        hours = seconds // 3600
        minutes = (seconds % 3600) // 60
        seconds = seconds % 60
        time_text = f"{days} days, {hours} hours, {minutes} minutes, {seconds} seconds"
        datatable = self.table()
        datatable.add_column("detail", justify="left", style="red", no_wrap=True)
        datatable.add_column("value", justify="left")
        datatable.add_row("Harpy:", f"v{__version__}")
        datatable.add_row("Time:", self.time_now())
        datatable.add_row("Duration:", time_text)
        datatable.add_row("Workflow Log: ", os.path.relpath(logfile))
        self.rule("[bold]Workflow Error[/]", style = "red")
        self.print(datatable)
        self.print("The workflow stopped due to an error. See the information Snakemake reported below.\n")


    def shell(self, text, rules: bool = False, style = None) -> None:
        """
        Prints the input text string as syntax-highlighted SHELL code to stderr
        """
        text = text.replace("(command exited with non-zero exit code)", "").strip()
        if rules:
            self.console.rule("Shell Code", style = 'dim')
        code = Syntax(text, "sh", background_color='default', dedent=True, code_width=2000, theme = "one-dark")
        self.print(code, soft_wrap=True, width = 2000, highlight = False, end = "", style = style)
        if rules:
            self.console.rule(style = 'dim')


    def validation(self, success: bool) -> None:
        '''
        If not `quiet`, print either a red x or green check following a validation log message
        '''
        if self.quiet == 0:
            if success:
                self.print("[green]🗸[/]")
            else:
                self.print("[red]𐄂[/]")


    def log(self, text, newline:bool = True):
        '''Print a rich-style log with the time in magenta and text in default'''
        _now =  time.strftime(r'[dim magenta]\[%H:%M:%S][/]')
        if self.quiet == 0:
            self.print(_now, text, highlight=False, end = "\n" if newline else " ")


    def progresspanel(self, progressbar: Progress, title: str|None = None, refresh: int = 2):
        """Returns a nicely formatted live-panel with the progress bar in it"""
        if self.quiet == 2:
            return nullcontext()
        return Live(
            Panel(
                progressbar, title = title, border_style="dim"
            ) if self.quiet != 2 else None,
            refresh_per_second=refresh,
            transient= self.quiet > 0,
            console=self.console
        )


    def progressbar(self) -> Progress:
        """
        The pre-configured transient progress bar that workflows and validations use
        """
        return Progress(
            TextColumn("{task.fields[active]}", style="yellow"),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(bar_width=None, complete_style="yellow", finished_style="dim blue"),
            TaskProgressColumn("{task.completed}/{task.total}", style = "blue") if self.quiet == 0 else TaskProgressColumn(style = "blue"),
            PausableTimeElapsedColumn(),
            transient = self.quiet > 0,
            auto_refresh = True,
            disable = self.quiet == 2,
            refresh_per_second=2,
            console= self.console,
            expand=True
        )


    def pulsebar(self, stderr: bool = False) -> Progress:
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
            disable = self.quiet == 2,
            console = self.console if stderr else None,
            expand=True
        )

    def process_sm_errors(self, errtext):
        '''
        final processing of the snakemake stderr text after an error has occured,
        returns early if ongoing or successful exit, otherwise processess the error text.
        If `isHPC=True`, will print the message: and reason: provided by snakemake, since
        it may point to an HPC-specific error not in the log file. 
        '''
        self.console.tab_size = 4
        self.console._highlight = False
        #print(errtext)
        #sys.exit()
        self.errortext = _Pushback(iter(errtext))
        self.missingoutput = []
        self.console.soft_wrap = True
        self.rules = []
        # shortcut to FileNotFoundError #
        line = next(self.errortext)
        if line.strip().startswith("FileNotFound"):
            if "envs/" in line and ".yaml'" in line:
                self.print("[red]Missing conda environment yaml file:[/][yellow]\n  " + line.split(":")[-1].replace("'", ""))
            else:
                self.print(line.strip(), style = "red")
            return
        # pick out conda-env errors
        if "Could not create conda" in line:
            for i in self.errortext:
                if "To search for alternate" in i:
                    break
                if i.strip():
                    if i.lstrip().startswith("-"):
                        self.print(i.rstrip(), soft_wrap = True, width = 2000, style = "bold red")
                    else:
                        self.print(i.rstrip(), soft_wrap = True, width = 2000, style = "red")
            return

        if ("Error" in line or "Exception" in line or "Missing input files" in line) and not ("RuleException" in line or "CalledProcessError" in line):
            #self.rule("[bold]Source of Error", style = "dim")
            self.print(line, highlight=False, soft_wrap = True, end = "", style = "red")
            for i in self.errortext:
                self.print(i, highlight = False, soft_wrap = True, end = "", style="red")
            return

        if "but some output files are missing" in line:
            self.missingoutput.append(line)
            for i in self.errortext:
                if "Shutting down, this might" in i:
                    break
                elif "but some output files are missing" not in i:
                    self.missingoutput[-1] += i
                else:
                    self.missingoutput.append(i)

        for i in self.errortext:
            if "Exiting because a job execution failed. Look below for error messages" in i:
                break
        for i in self.errortext:
            if "(100%) done" in i:
                break
            if "Error in group" in i:
                self.rule("[bold]Source of Error", style="dim")
                i = next(self.errortext).strip()
            if i.startswith("[") and i.strip().endswith("]"):
                break

        # error in rule line — now populates SnakeRule instead of printing
        for i in self.errortext:
            if "(100%) done" in i:
                break
            if "RuleException" in i:
                sys.exit(1)
            m = _HEADER_RE.match(i.strip())
            if m:
                rule = self._parse_rule_block()
                rule.name = m.group(1)
                self.rules.append(rule)
            elif i.startswith("Complete log"):
                #return
                break
            elif i.startswith("WorkflowError"):
                break #return

        if self.missingoutput:
            self.print("\n[bold dim]──── ⚠ Error Reported by Snakemake")
            self.print(format_missing_output_block(self.missingoutput[0]), style = "red")
        else:
            self.print_ruleerror(rule)

    def _parse_rule_block(self) -> SnakeRule:
        rule = SnakeRule()
        key = None
        for line in self.errortext:
            if not line.strip():
                continue
            m = _KEY_RE.match(line)
            if m:
                key, val = m.group(1), m.group(2)
                if key == "jobid":
                    rule.jobid = int(val)
                elif key == "input":
                    rule.input = val
                elif key == "output":
                    rule.output = val
                elif key == "log":
                    pass
                    #rule.logfile = rule.logfile or val
                    #skip in favor of parsing logs later
                    #rule.logs.append(val)
                elif key == "message":
                    rule.message = val
                elif key == "conda-env":
                    rule.env = val
                elif key == "resources":
                    rule.resources = [r.strip() for r in val.split(",")]
                elif key == "shell":
                    rule.cmd = ""
            elif key == "shell" and line[:1].isspace():
                rule.cmd += line + "\n"
            elif line[:1].isspace():
                continue
            else:
                self.errortext.push(line)
                break
        self._consume_logfile_blocks(rule)
        rule.cmd = format_shell_cmd(rule.cmd)
        return rule

    def _consume_logfile_blocks(self, rule: SnakeRule):
        """Snakemake prints 'Logfile <path>:' + fenced content, or
        'Logfile <path> ... not found.', right after a rule's error block.
        Consumes zero or more of these, pushing back the first line that isn't one."""
        for line in self.errortext:
            m = _LOGFILE_HEADER_RE.match(line.strip())
            if m:
                content = []
                for j in self.errortext:
                    if _is_separator(j):
                        if content:
                            break  # closing fence
                        continue  # opening fence
                    content.append(j)

                logtext = escape(re.sub(r'\n{3,}', '\n\n', "".join(content)).removeprefix("    "))
                # purge out all unnecessary papermill error text
                if "papermill.exceptions.PapermillExecutionError:" in logtext:
                    logtext = logtext.partition("papermill.exceptions.PapermillExecutionError:")[-1]
                    chunks = logtext.split("\n\n")
                    filtered = [c for c in chunks if not c.startswith("File ")]
                    logtext = "\n\n".join(filtered)
                rule.logs[m.group(1)] = logtext.strip()
                continue
            if _LOGFILE_NOTFOUND_RE.match(line.strip()):
                continue  # path already in rule.logs, nothing to capture
            self.errortext.push(line)
            break

    def print_ruleerror(self, rule: SnakeRule):
        self.rule(f"[default bold]Triggering Rule[/][yellow bold] {rule.name}", style = "yellow")

        self.print("input:", style = "bold default")
        self.print("  " + rule.input.replace(", ", "\n  "), style = "red")

        self.print("output:", style = "bold default")
        self.print("  " + rule.output.replace(", ", "\n  "), style = "red")

        if rule.logs:
            self.print("log:", style = "bold default")
            for i in rule.logs:
                _i = i.replace("(check log file(s) for error details)", "").strip()
                self.print("  " + _i, style = "red")

        if rule.env:
            env_clean = rule.env.split("/")[:-2]
            self.print(
                f"[bold default]conda-env:[/] [red]{'/'.join(env_clean)}[/]"
            )

        if rule.message != "None":
            self.print("Message:", style = "bold default")
            self.print(rule.message, style = "red")

        if rule.cmd:
            self.print("")
            self.print("[bold dim]──── ❯ Command Invoked")
            self.shell(rule.cmd)

        if rule.logs:
            self.print("")

        for name, content in rule.logs.items():
            self.print(f"──── 🗎 {name}", style = "bold dim")
            self.print(content, style = "red")

    def format_shell(self):
        '''format the snakemake rule shell command nicely and print it to the console'''
        text = ""
        for i in self.errortext:
            if "(command exited" in i:
                break
            text += i
        self.print("")
        self.print("[bold dim]──── ❯ Command Invoked")
        self.shell(text.strip("\n"))


    def print_logfile(self, errline):
        '''process and print the contents of a logfile in the snakemake error log'''
        merged_text = ""
        _log = errline.rstrip().split()[1]
        self.print("")
        self.print(f"[bold dim]──── 🗎 {_log.rstrip(':')}")
        if "empty file" in errline:
            self.print(f"{_log.replace(':','')} is empty\n", style = "dim")
            _ = next(self.errortext)
            return
        if "not found" in errline:
            self.print(f"{_log} was not found\n", style = "red")
            _ = next(self.errortext)
            return
        lines = 0
        for i in self.errortext:
            if lines == 2:
                break
            if "====" in i:
                lines += 1
                continue
            merged_text += i
        logtext = escape(re.sub(r'\n{3,}', '\n\n', merged_text).removeprefix("    "))
        # purge out all unnecessary papermill error text
        if "papermill.exceptions.PapermillExecutionError:" in logtext:
            logtext = logtext.partition("papermill.exceptions.PapermillExecutionError:")[-1]
            chunks = logtext.split("\n\n")
            filtered = [c for c in chunks if not c.startswith("File ")]
            logtext = "\n\n".join(filtered)
        self.print("[red]" + logtext, overflow = "ignore", crop = False)
        _ = next(self.errortext)
