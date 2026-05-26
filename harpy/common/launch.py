"""launch snakemake"""

from datetime import datetime
import os
import re
import signal
import sys
import subprocess
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
        """
        Initialize the launcher, start Snakemake, and run setup, startup checks, and job monitoring.
        
        This constructor starts the Snakemake subprocess with the provided arguments, initializes internal tracking state (timing, exit codes, job inventory, progress UI, and error collection), installs foreground/background signal handlers, and drives the launch sequence by calling workflow_setup(), check_startup(), and monitor_jobs() in order. It handles user interrupt (prints a termination message and exits with status 1), and on completion or error it stops the progress UI, removes signal handlers, collects remaining stderr into the error log, ensures the subprocess is terminated, and purges empty log files in the output directory.
        
        Parameters:
            sm_args (str): Command-line arguments to launch Snakemake (space-separated).
            outdir (str): Path to the output directory where logs may be created.
            quiet (int): Verbosity level (higher values suppress more interactive output).
            printer (HarpyPrint): Console/printer instance used for progress and messages.
        """
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
                    self.process.communicate(timeout=5)
                except subprocess.TimeoutExpired:
                    self.process.kill()
                    self.process.communicate()
            purge_empty_logs(outdir)

    def _is_foreground(self) -> bool:
        """
        Determine whether the current process is in the terminal's foreground process group.
        
        Returns:
            bool: `True` if the current process group equals the terminal's foreground process group, `False` otherwise (returns `False` if the foreground group cannot be queried).
        """
        try:
            return os.getpgrp() == os.tcgetpgrp(sys.stdin.fileno())
        except OSError:
            return False

    def _handle_sigtstp(self, signum, frame):
        """
        Handle SIGTSTP (Ctrl+Z) by silencing live console output and then suspending the process.
        
        Flushes any pending file output and sets the console into quiet mode so live progress output stops. Restores the SIGTSTP handler to the default and re-raises SIGTSTP so the operating system actually suspends the process.
        
        Parameters:
            signum (int): Signal number (SIGTSTP).
            frame (types.FrameType): Current stack frame at the time the signal was received.
        """
        self.print.file.flush()
        self.print.console.quiet = True
        # reset to default so the process actually suspends
        signal.signal(signal.SIGTSTP, signal.SIG_DFL)
        signal.raise_signal(signal.SIGTSTP)

    def _handle_sigcont(self, signum, frame):
        """
        Restore live console output when the process resumes in the terminal foreground.
        
        Re-enables the console (clears quiet mode) and emits a blank line to reset the cursor if the current
        process group is the terminal foreground. Also reinstates the SIGTSTP handler.
        
        Parameters:
            signum (int): Signal number received.
            frame (frame): Current stack frame at signal delivery.
        """
        signal.signal(signal.SIGTSTP, self._handle_sigtstp)
        if self._is_foreground():
            self.print.console.quiet = False
            self.print.console.print("")  # force cursor to a fresh line

    def _setup_bg_signal_handlers(self):
        """
        Install background/foreground signal handlers to manage terminal suspend/resume behavior.
        
        Sets handlers for SIGTSTP and SIGCONT so the launcher can quiet console output when suspended and restore it when resumed.
        """
        signal.signal(signal.SIGTSTP, self._handle_sigtstp)
        signal.signal(signal.SIGCONT, self._handle_sigcont)

    def _teardown_bg_signal_handlers(self):
        """
        Restore the default handling for terminal stop (SIGTSTP) and continue (SIGCONT) signals.
        
        This returns the process's SIGTSTP and SIGCONT handlers to the system defaults, undoing any custom handlers previously installed.
        """
        signal.signal(signal.SIGTSTP, signal.SIG_DFL)
        signal.signal(signal.SIGCONT, signal.SIG_DFL)

    def is_done(self) -> bool:
        """
        Determine whether the Snakemake subprocess has finished or an exit code has been recorded.
        
        Returns:
            true if an exit code has been set (exitcode > -1) or the subprocess has terminated, false otherwise.
        """
        if self.exitcode > -1 or self.process.poll():
            return True
        return False

    def update_total_active(self):
        '''update self.total_active with the sum of all the active jobs in self.rule_inventory'''
        self.total_active = sum(self.job_inventory[rule].active() for rule in self.job_inventory if rule != "total")

    def nothing_to_do(self):
        """
        Terminate the launcher early when Snakemake reports there is no work to perform.
        
        If the current captured stderr line contains the substring "Nothing to be", prints
        "All outputs already present" in green and exits the process with status 0.
        
        Raises:
            SystemExit: exits with status code 0 when no work is detected.
        """
        if "Nothing to be" in self.output:
            self.print.rule("[bold]All outputs already present", style="green")
            sys.exit(0)

    def iserror(self) -> bool:
        """
        Determine if the current captured stderr line contains Snakemake error-trigger words.
        
        Returns:
            `True` if the current output contains any of the strings "Exception", "Error", or "MissingOutputException", `False` otherwise.
        """
        return "Exception" in self.output or "Error" in self.output or "MissingOutputException" in self.output

    def nextline(self, strip: bool = False):
        """reads the next line of stderr"""
        _ = self.process.stderr.readline()
        if not _:
            self.output = ""
        else:
            self.output = _.strip() if strip else _

    def pause_progress(self, rulename):
        """
        Pause the elapsed-time column for a rule's progress task.
        
        Parameters:
            rulename (str): The rule name whose progress task's elapsed-time column will be paused.
        """
        self.progress.columns[4].pause(self.task_ids[rulename])

    def resume_progress(self, rulename):
        """
        Resume the elapsed-time column for a rule's progress task.
        
        Parameters:
            rulename (str): The rule identifier whose progress timer should be resumed.
        """
        self.progress.columns[4].resume(self.task_ids[rulename])

    def update_finished_progress(self):
        """
        Update progress tasks when a job completion line appears in stderr.
        
        Extract the completed job id from self.output, remove that id from the corresponding rule's tracked ids, recalculate total active jobs, advance the per-rule and aggregate progress tasks, pause the rule's elapsed-time column when it has no remaining active jobs, and mark tasks as finished or clear the total-active label when their totals complete.
        """
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
        """
        Check Snakemake startup stderr for early termination or error conditions and update the launch exit state.
        
        Reads the next stderr line; if the subprocess has already exited or the line indicates an error, set self.exitcode to success when the process exited cleanly or to the snakefile error code otherwise. If the output contains "Conda", override the exit code to the conda-related error. After setting the code, drain and print any remaining stderr lines to the console in red.
        """
        self.nextline()
        if self.process.poll() or self.iserror():
            self.exitcode = EXIT_CODE_SUCCESS if self.process.poll() == 0 else EXIT_CODE_SNAKEFILE_ERROR
            self.exitcode = EXIT_CODE_CONDA_ERROR if "Conda" in self.output else self.exitcode
            while self.output:
                self.print.print(self.output, style="red")
                self.nextline()

    def workflow_setup(self):
        """
        Process Snakemake's startup stderr until the job summary table is reached and prepare launcher state for runtime monitoring.
        
        This consumes stderr lines produced during workflow setup, detects early terminal conditions (e.g., "Nothing to be", missing inputs, errors, or dependency installation steps), and updates launcher state accordingly. Side effects include:
        - setting self.exitcode for startup or snakemake errors,
        - calling sys.exit(0) when no work is needed,
        - appending startup error lines to self.errorlog,
        - setting self.deps and augmenting self.deploy_text when dependency/build steps are detected,
        - populating self.job_inventory with Rule entries parsed from the "Job stats:" table and adjusting the aggregate total.
        """
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
        """
        Monitor Snakemake stderr and update progress UI and internal state during workflow execution.
        
        Reads Snakemake stderr lines until completion or error, detecting rule start lines to create or update per-rule progress tasks, recording job IDs in self.job_inventory[rule].ids, and detecting finished job lines to advance progress. Sets self.exitcode on runtime errors or when execution finishes, and updates the aggregate total progress task in self.task_ids.
        """
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
                if self.iserror() or self.process.poll() == 1:
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
        """
        Normalize a non-negative exit code or collect remaining stderr lines into the instance error log.
        
        If self.exitcode is less than or equal to 0, set it to 0 and return. Otherwise, read the remainder of the subprocess stderr stream and append each line to self.errorlog except lines that end with ", in <module>".
        """
        if self.exitcode <= 0:
            self.exitcode = max(self.exitcode, 0)
            return
        stderr_output = self.process.stderr.read()
        for line in stderr_output.splitlines(keepends=True):
            if not line.strip().endswith(", in <module>"):
                self.errorlog.append(line)