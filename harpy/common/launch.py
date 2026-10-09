"""launch snakemake"""

import os
import re
import signal
import subprocess
import sys
from datetime import datetime

from harpy.common.file_ops import purge_empty_logs
from harpy.common.printing import HarpyPrint
from harpy.common.progress import PanelProgress

EXIT_CODE_SUCCESS = 0
EXIT_CODE_SNAKEFILE_ERROR = 1
EXIT_CODE_CONDA_ERROR = 2
EXIT_CODE_RUNTIME_ERROR = 3
# quiet = 0 : print all things, full progressbar
# quiet = 1 : print all text, only "Total" progressbar
# quiet = 2 : print nothing, no progressbar

# snakemake logs the error of a job that has retries left, followed immediately by "Trying to restart job N."
_RESTART_RE = re.compile(r"^Trying to restart job (\d+)\.")
_ERROR_HEADER_RE = re.compile(r"^Error in (?:rule|group) ")
_JOBID_RE = re.compile(r"^\s+jobid:\s*(\d+)")
_ALL_DONE_RE = re.compile(r"^\d+ of \d+ steps \(100%\) done")
# what snakemake prints once it has decided that the workflow failed
_FAILED_MARKERS = ("Shutting down, this might take some time.", "Exiting because a job execution failed")

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
        # a job error whose outcome isn't known yet: snakemake retries jobs that have retries left (--retries)
        self.pending_error: list[str] | None = None
        self.error_line: str = ""
        self.pending_jobs: set[int] = set()
        self.pending_attributed: bool = False
        self.awaiting_jobid: bool = False
        self.retries: dict[int, int] = {}
        self.output: str = ""
        self.job_inventory: dict = {}
        self.task_ids: dict = {}
        self.total_active: int = 0
        self.print = printer
        self._setup_bg_signal_handlers()
        self.progress = PanelProgress(console=self.print.console, quiet=quiet, transient = quiet==2 ).bar()
        #self.progress = self.print.progressbar()

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

    def job_failed(self) -> bool:
        '''
        Decide whether the line in `self.output` means the workflow has failed, tolerating retries.

        Snakemake prints the error of a failed job, and then either "Trying to restart job N." if it has retries
        left, or later "Shutting down..." when it gives up. So the first line that looks like an error only starts a
        *pending* error: its text is kept, and the verdict comes from what follows (progress lines keep being
        processed in the meantime, since other jobs go on and can even finish between the error and the verdict).
        - a restart of every job with a pending error means there is nothing wrong (yet), the pending error is dropped
        - the failure markers, or the process exiting with an error, mean it has failed. In that case `self.output`
          and `self.errorlog` are set as if the first error line had stopped monitoring, which is what the error
          parsing expects.
        Assumes that the output of the jobs themselves doesn't end up in the output of snakemake, which is the case for
        the scheduler plugins (their jobs write to log files) but not for e.g. a blocking submit command that passes it on.
        '''
        line = self.output
        if self.pending_error is None:
            if self.iserror():
                self.pending_error = []
                self.error_line = line
                self.pending_jobs = set()
                self.pending_attributed = bool(_ERROR_HEADER_RE.match(line.strip()))
                self.awaiting_jobid = self.pending_attributed
                return self._exited_with_error(line) and self._confirm_failure()
            return self._exited_with_error(line)

        self.pending_error.append(line)
        if _ERROR_HEADER_RE.match(line.strip()):
            self.pending_attributed = True
            self.awaiting_jobid = True
        elif self.awaiting_jobid and (m := _JOBID_RE.match(line)):
            self.pending_jobs.add(int(m.group(1)))
            self.awaiting_jobid = False
        elif (m := _RESTART_RE.match(line)):
            jobid = int(m.group(1))
            self.retries[jobid] = self.retries.get(jobid, 0) + 1
            self.pending_jobs.discard(jobid)
            if self.pending_attributed and not self.pending_jobs:
                # every failed job is being retried, forget about it
                self.pending_error = None
                self.error_line = ""
                return False
        if line.startswith(_FAILED_MARKERS) or self._exited_with_error(line):
            return self._confirm_failure()
        return False

    def _exited_with_error(self, line: str) -> bool:
        '''
        Whether snakemake has exited with an error, and all of its output has been read. The exit code alone isn't a verdict:
        snakemake can exit while there is still output in the pipe, e.g. when this is slower than snakemake is, and what is
        left can decide the outcome (a job being retried, the last progress updates).
        '''
        # at the end of the output the process is done or about to be, so wait for its exit code rather than poll (which can still say None)
        return not line and self.process.wait() != 0

    def _confirm_failure(self) -> bool:
        '''A pending error turned out to be a failure: hand its text over as if the first error line had stopped monitoring'''
        if self.pending_error is not None:
            self.output = self.error_line
            self.errorlog.extend(i for i in self.pending_error if not i.strip().endswith(", in <module>"))
            self.pending_error = None
        return True

    def eof_failed(self) -> bool:
        '''At the end of the output, whether a pending error was a failure. Settled by the exit code of snakemake.'''
        if self.process.wait() != 0:
            return self._confirm_failure()
        # exited without an error, so it was only text that looked like one
        self.pending_error = None
        return False

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
        if self.quiet == 2:
            return
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
                #progress = PanelProgress(self.print.console, self.quiet, title=self.deploy_text).pulse()
                with PanelProgress(self.print.console, self.quiet, title=self.deploy_text, transient=True).pulse() as progress:
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
                    #progress.stop()
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
        with self.progress:
            self.task_ids["total_progress"] = self.progress.add_task(
                "[bold blue]Progress",
                total=self.job_inventory["total"].total,
                active="[bold]0"
            )
            while self.output:
                self.nextline()
                if self.job_failed():
                    self.exitcode = EXIT_CODE_RUNTIME_ERROR
                    break
                # while an error is pending, its text (error blocks, log contents) is also going through here, so be strict
                pending = self.pending_error is not None
                if (_ALL_DONE_RE.match(self.output) if pending else "(100%) done" in self.output) or self.output.startswith("Nothing to be"):
                    self.exitcode = EXIT_CODE_SUCCESS
                    break
                if self.output.startswith("Complete log") or self._exited_with_error(self.output):
                    self.exitcode = EXIT_CODE_SUCCESS if self.process.poll() == 0 else EXIT_CODE_RUNTIME_ERROR
                    break
                # (group error blocks list their rules, indented, in the same way group jobs are started)
                if (self.output if pending else self.output.lstrip()).startswith(("rule ", "localrule ")):
                    rule = self.output.split()[-1].replace(":", "")
                    if rule not in self.task_ids and rule != "all":
                        self.task_ids[rule] = self.progress.add_task(self.job_inventory[rule].name, total=self.job_inventory[rule].total, visible=self.quiet != 1, active=1)
                    while True:
                        self.nextline()
                        if not self.output:                              # EOF: process died
                            if self.pending_error is None or self.eof_failed():
                                self.exitcode = EXIT_CODE_RUNTIME_ERROR
                            return
                        if self.job_failed():
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
            # out of output without a verdict: only a pending error is left to settle, by how snakemake exited
            if self.exitcode < 0 and self.pending_error is not None and self.eof_failed():
                self.exitcode = EXIT_CODE_RUNTIME_ERROR

    def return_or_collect(self):
        if self.exitcode <= 0:
            self.exitcode = max(self.exitcode, 0)
            return
        for line in self.process.stderr:
            if not line.strip().endswith(", in <module>"):
                self.errorlog.append(line)