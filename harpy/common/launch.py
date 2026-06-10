"""launch snakemake"""

import os
import re
import signal
import subprocess
import sys
from datetime import datetime

from harpy.common.file_ops import purge_empty_logs
from harpy.common.printing import HarpyPrint

EXIT_CODE_SUCCESS = 0
EXIT_CODE_SNAKEFILE_ERROR = 1
EXIT_CODE_CONDA_ERROR = 2
EXIT_CODE_RUNTIME_ERROR = 3
# quiet = 0 : print all things, full progressbar
# quiet = 1 : print all text, only "Total" progressbar
# quiet = 2 : print nothing, no progressbar

class Rule:
    """A class that stores job information with the fields: name, total, ids"""
    def __init__(self, name, total):
        self.name: str = name
        self.total: int = total
        self.ids: set = set()

    def active(self) -> int:
        return len(self.ids)

class LaunchSnakemake():
    """launch snakemake with the given commands and monitor its progress"""
    def __init__(self, sm_args, outdir, quiet, printer: HarpyPrint):
        self.exitcode = -1
        self.start_time = datetime.now()
        self.deps: bool = False
        self.deploy_text: str = ""
        self.quiet = quiet
        self.cmd: list[str] = sm_args.split()
        self.outdir: str = outdir
        self.process = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.errorlog = []
        self.output: str = ""
        self.job_inventory: dict = {}
        self.task_ids: dict = {}
        self.total_active: int = 0
        self.print = printer
        self._setup_bg_signal_handlers()
        self.progress = self.print.progressbar()

        try:
            self.workflow_setup()
            if self.is_done():
                return
            self.check_startup()
            if self.is_done():
                return
            self.monitor_jobs()
        except KeyboardInterrupt:
            if self.quiet < 2:
                for _ in range(1):
                    self.print.file.write("\033[F\033[K")
                self.print.file.flush()
            self.print.print("")
            self.print.rule("[bold]Terminating Harpy", style="yellow")
            sys.exit(1)
        finally:
            self.progress.stop()
            self._teardown_bg_signal_handlers()
            self.return_or_collect()
            if self.process.poll() is None:
                self.process.terminate()
                try:
                    self.process.communicate(timeout=0.5)
                except subprocess.TimeoutExpired:
                    self.process.kill()
                    #self.process.communicate()
            purge_empty_logs(outdir)

    def _is_foreground(self) -> bool:
        """Check if this process is in the terminal's foreground process group."""
        try:
            return os.getpgrp() == os.tcgetpgrp(sys.stdin.fileno())
        except OSError:
            return False

    def _handle_sigtstp(self, signum, frame):
        """Suppress Live output when backgrounded via Ctrl+Z."""
        self.print.file.flush()
        self.print.console.quiet = True
        # reset to default so the process actually suspends
        signal.signal(signal.SIGTSTP, signal.SIG_DFL)
        signal.raise_signal(signal.SIGTSTP)

    def _handle_sigcont(self, signum, frame):
        """Restore Live output when foregrounded again, but only if in foreground."""
        signal.signal(signal.SIGTSTP, self._handle_sigtstp)
        if self._is_foreground():
            self.print.console.quiet = False
            self.print.console.print("")  # force cursor to a fresh line

    def _setup_bg_signal_handlers(self):
        '''Setup a signal handler to make sure the progress bar vanishes if harpy is put in the backgroud'''
        signal.signal(signal.SIGTSTP, self._handle_sigtstp)
        signal.signal(signal.SIGCONT, self._handle_sigcont)

    def _teardown_bg_signal_handlers(self):
        '''Restore signal handlers'''
        signal.signal(signal.SIGTSTP, signal.SIG_DFL)
        signal.signal(signal.SIGCONT, signal.SIG_DFL)

    def is_done(self) -> bool:
        '''check if self.exitcode > -1 or a value exists for self.process.poll()'''
        if self.exitcode > -1 or self.process.poll():
            return True
        return False

    def update_total_active(self):
        '''update self.total_active with the sum of all the active jobs in self.rule_inventory'''
        self.total_active = sum(self.job_inventory[rule].active() for rule in self.job_inventory if rule != "total")

    def nothing_to_do(self):
        '''check if self.output has the "Nothing to be" triggering text and exit with return code 0 if true'''
        if "Nothing to be" in self.output:
            self.print.rule("[bold]All outputs already present", style="green")
            sys.exit(0)

    def iserror(self) -> bool:
        '''logical check for erroring trigger words in snakemake output'''
        return "Exception" in self.output or "Error" in self.output or "MissingOutputException" in self.output

    def nextline(self, strip: bool = False):
        """reads the next line of stderr"""
        _ = self.process.stderr.readline()
        if not _:
            self.output = ""
        else:
            self.output = _.strip() if strip else _

    def pause_progress(self, rulename):
        '''pause the time elapsed col for a rule's progress bar'''
        self.progress.columns[4].pause(self.task_ids[rulename])

    def resume_progress(self, rulename):
        '''resume the time elapsed col for a rule's progress bar'''
        self.progress.columns[4].resume(self.task_ids[rulename])

    def update_finished_progress(self):
        '''Process the stderr output and update the progressbars accordingly'''
        completed = int(re.search(r"\d+", self.output).group())
        for job, details in self.job_inventory.items():
            if completed in details.ids:
                self.job_inventory[job].ids.discard(completed)
                self.update_total_active()
                task_id = self.task_ids[job]
                _active = self.job_inventory[job].active()
                if _active < 1:
                    self.pause_progress(job)
                    self.progress.update(task_id, advance=1, refresh=True, active="[dim yellow]⋯")
                else:
                    self.progress.update(task_id, advance=1, refresh=True, active=_active)
                self.progress.update(self.task_ids["total_progress"], refresh=True, advance=1, active=f"[bold]{self.total_active}")
                if self.progress.tasks[self.task_ids[job]].completed == self.progress.tasks[task_id].total:
                    self.progress.update(self.task_ids[job], refresh=True, description=f"[dim]{details.name}", active="[dim blue]✓")
                break
        if self.progress.tasks[self.task_ids["total_progress"]].completed == self.progress.tasks[self.task_ids["total_progress"]].total:
            self.progress.update(self.task_ids["total_progress"], refresh=True, active=" ")

    def check_startup(self):
        '''monitors the process for startup errors or things already being done'''
        self.nextline()
        if self.process.poll() or self.iserror():
            self.exitcode = EXIT_CODE_SUCCESS if self.process.poll() == 0 else EXIT_CODE_SNAKEFILE_ERROR
            self.exitcode = EXIT_CODE_CONDA_ERROR if "Conda" in self.output else self.exitcode
            while self.output:
                self.print.print(self.output, style="red")
                self.nextline()

    def workflow_setup(self):
        '''processes the workflow setup text snakemake prints to the console up to the end of the job summary table'''
        while self.exitcode < 0:
            if self.quiet < 2:
                with self.print.status("[dim]Preparing workflow", spinner="point", spinner_style="yellow"):
                    while self.output.startswith("Building DAG of jobs...") or self.output.startswith("Assuming"):
                        self.nextline()
                if "Nothing to be" in self.output:
                    self.print.rule("[bold]All outputs already present", style="green")
                    sys.exit(0)
            else:
                while self.output.startswith("Building DAG of jobs...") or self.output.startswith("Assuming"):
                    self.nextline()
            while not self.output.startswith("Job stats:") and self.exitcode < 0:
                if "Creating conda environment" in self.output or "Running post-deploy" in self.output:
                    self.deps = True
                    self.deploy_text += "[dim]Installing workflow software"
                    break
                if "Pulling singularity image" in self.output:
                    self.deps = True
                    self.deploy_text += "[dim]Building software container"
                    break
                if "Nothing to be" in self.output:
                    self.print.rule("[bold]All workflow outputs already present", style="green")
                    sys.exit(0)
                if "MissingInput" in self.output:
                    self.exitcode = EXIT_CODE_SNAKEFILE_ERROR
                    return
                if "Error" in self.output or "Exception" in self.output:
                    self.exitcode = EXIT_CODE_SNAKEFILE_ERROR
                    self.errorlog.append(self.output)
                    return
                self.nextline()
            if self.deps:
                progress = self.print.pulsebar()
                with self.print.progresspanel(progress, title=self.deploy_text, refresh=8):
                    _taskid = progress.add_task("[dim]Working...", total=None)
                    while not self.output.startswith("Job stats:"):
                        if "Creating conda environment" in self.output:
                            _desc = self.output.split()[-1].removesuffix("...")
                            progress.update(_taskid, description=_desc)
                        self.nextline()
                        if self.process.poll() or self.iserror():
                            self.exitcode = EXIT_CODE_SUCCESS if self.process.poll() == 0 else 2
                            break
                        self.nothing_to_do()
                    progress.stop()
            if self.process.poll() or self.exitcode >= 0:
                return
            self.nothing_to_do()
            while True:
                self.nextline()
                if self.output.startswith("Select jobs to execute"):
                    self.job_inventory["total"].total -= 1
                    return
                try:
                    rule, count = self.output.split()
                    if rule in ["job", "all"] or "----" in rule:
                        continue
                    rule_desc = rule.replace("_", " ")
                    self.job_inventory[rule] = Rule(rule_desc, int(count))
                except ValueError:
                    pass
                if self.process.poll() or self.iserror():
                    self.exitcode = EXIT_CODE_SUCCESS if self.process.poll() == 0 else EXIT_CODE_SNAKEFILE_ERROR
                    return

    def monitor_jobs(self):
        '''monitors the Snakemake stderr output while jobs are running'''
        if self.is_done():
            return
        with self.print.progresspanel(self.progress):
            self.task_ids["total_progress"] = self.progress.add_task(
                "[bold blue]Progress",
                total=self.job_inventory["total"].total,
                active="[bold]0"
            )
            while self.output:
                self.nextline()
                if self.iserror() or (self.process.poll() not in (None, 0)):
                    self.exitcode = EXIT_CODE_RUNTIME_ERROR
                    break
                if "(100%) done" in self.output or self.output.startswith("Nothing to be") or self.process.poll() == 0:
                    self.exitcode = EXIT_CODE_SUCCESS
                    break
                if self.output.startswith("Complete log") or self.process.poll():
                    self.exitcode = EXIT_CODE_SUCCESS if self.process.poll() == 0 else EXIT_CODE_RUNTIME_ERROR
                    break
                if self.output.lstrip().startswith("rule ") or self.output.lstrip().startswith("localrule "):
                    rule = self.output.split()[-1].replace(":", "")
                    if rule not in self.task_ids and rule != "all":
                        self.task_ids[rule] = self.progress.add_task(self.job_inventory[rule].name, total=self.job_inventory[rule].total, visible=self.quiet != 1, active=1)
                    while True:
                        self.nextline()
                        if not self.output:                              # EOF: process died
                            self.exitcode = EXIT_CODE_RUNTIME_ERROR
                            return
                        if self.iserror() or self.process.poll() not in (None, 0):
                            self.exitcode = EXIT_CODE_RUNTIME_ERROR
                            return
                        if "jobid: " in self.output:
                            job_id = int(self.output.strip().split()[-1])
                            if rule != "all":
                                self.job_inventory[rule].ids.add(job_id)
                                self.resume_progress(rule)
                                self.progress.update(self.task_ids[rule], active=self.job_inventory[rule].active())
                                self.update_total_active()
                                self.progress.update(self.task_ids["total_progress"], refresh=True, active=f"[bold]{self.total_active}")
                            break
                if self.output.startswith("Finished jobid: "):
                    self.update_finished_progress()

    def return_or_collect(self):
        if self.exitcode <= 0:
            self.exitcode = max(self.exitcode, 0)
            return
        for line in self.process.stderr:
            if not line.strip().endswith(", in <module>"):
                self.errorlog.append(line)
