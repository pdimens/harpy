"""launch snakemake"""

from datetime import datetime
import re
import os
import sys
import subprocess
from rich import box, console
from rich.syntax import Syntax
from rich.table import Table
from harpy.common.file_ops import gzip_file, purge_empty_logs
from harpy.common.printing import CONSOLE, print_onerror, print_setup_error
from harpy.common.progress import harpy_progressbar, harpy_pulsebar, harpy_progresspanel

EXIT_CODE_SUCCESS = 0
EXIT_CODE_GENERIC_ERROR = 1
EXIT_CODE_CONDA_ERROR = 2
EXIT_CODE_RUNTIME_ERROR = 3
# quiet = 0 : print all things, full progressbar
# quiet = 1 : print all text, only "Total" progressbar
# quiet = 2 : print nothing, no progressbar

class launch_snakemake():
    """launch snakemake with the given commands"""
    def __init__(self, sm_args, outdir, sm_logfile, quiet, CONSOLE = CONSOLE):
        self.exitcode = -1
        self.start_time = datetime.now()
        self.deps: bool = False
        self.deploy_text: str = ""
        self.logfile = sm_logfile
        self.quiet = quiet
        self.cmd: list[str] = sm_args.split()
        self.outdir: str = outdir
        self.process = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
        self.output: str = ""
        self.job_inventory: dict = {}
        self.total_active: int = 0

        try:
            self.workflow_setup()
            if self.exitcode > -1 or self.process.poll():
                self.process_finish()
            self.check_startup()
            if self.exitcode > -1 or self.process.poll():
                if self.exitcode == 0:
                    return
                else:
                    self.process_finish()
                    return
            self.monitor_jobs()
            #CONSOLE.print(self.exitcode)
            #if self.exitcode > -1 or self.process.poll():
                #CONSOLE.print("AFTERMONITOR")
            self.process_finish()
        except KeyboardInterrupt:
            CONSOLE.print("")
            CONSOLE.rule("[bold]Terminating Harpy", style = "yellow")
            sys.exit(1)
        finally:
            self.process.terminate()
            self.process.wait()
            if os.path.exists(os.path.join(outdir,sm_logfile)):
                gzip_file(os.path.join(outdir,sm_logfile))
            purge_empty_logs(outdir)

    def update_total_active(self):
        self.total_active = sum([len(self.job_inventory[rule][2]) for rule in self.job_inventory if rule != "total"])

    def iserror(self) -> bool:
        """logical check for erroring trigger words in snakemake output"""
        return "Exception" in self.output or "Error" in self.output or "MissingOutputException" in self.output

    def process_error(self):
        '''interpret rule errors and print them with nice format'''
        if self.output.strip().startswith("Logfile"):
            self.print_logfile()
            return
        text = self.output.removeprefix("    ").rstrip().lstrip()
        valid_keys = ["jobid","input","output","log","conda-env","container","shell","wildcards", "affected files"]
        _split = text.split(':')
        if _split[0] == "message":
            return
        if len(_split) == 1 or _split[0] not in valid_keys:
            CONSOLE.print(f"[red]{text}", overflow = "ignore", crop = False)
            return
        key = _split[0]
        vals = [i.strip() for i in _split[1].split(",")]
        if len(vals) == 1:
            CONSOLE.print(f"[bold default]{key}: [/][red]" + "".join(vals))
            return
        CONSOLE.print(f"[bold default]{key}: [/]\n  [red]" + "\n  ".join(vals))

    def print_shellcmd(self):
        _table = Table(
            show_header=False,
            pad_edge=False,
            show_edge=False,
            padding=(0,0),
            box=box.SIMPLE,
        )
        _table.add_column("Lpadding", justify="left")
        _table.add_column("shell", justify="left")
        _table.add_column("Rpadding", justify="left")

        text = ""
        while "(command exited" not in self.output or not self.output:
            self.nextline()
            if not self.output:
                break
            text += self.output

        text = text.replace("(command exited with non-zero exit code)", "").rstrip().lstrip().replace("\t", "  ")
        text = re.sub(r' {2,}|\t+', '  ', text)
        cmd = Syntax(text, lexer = "bash", tab_size=2, word_wrap=True, padding=1, dedent=True, theme = "paraiso-dark")
        _table.add_row("  ", cmd, "  ")
        CONSOLE.print("[bold default]shell:", _table)

    def print_logfile(self):
        '''process and print the contents of a logfile in the snakemake error log'''
        merged_text = ""
        _log = self.output.rstrip().split()[1]
        CONSOLE.rule(f"[bold]Log File: {_log.rstrip(':')}", style = "yellow")
        lines = 0
        while lines < 2:
            self.nextline()
            if "====" in self.output:
                lines += 1
                continue
            merged_text += self.output
        if "====" in self.output:
            CONSOLE.print("[red]" + re.sub(r'\n{3,}', '\n\n', merged_text), overflow = "ignore", crop = False)
            self.nextline()

    def nextline(self, strip: bool = False):
        """reads the next line of stderr"""
        _ = self.process.stderr.readline()
        if not _:
            self.output = ""
        else:
            self.output = _.strip() if strip else _

    def check_startup(self):
        """monitors the process for startup errors or things already being done"""
        self.nextline()
        # check for syntax errors at the very beginning
        if self.process.poll() or self.iserror():
            self.exitcode = EXIT_CODE_SUCCESS if self.process.poll() == 0 else EXIT_CODE_GENERIC_ERROR
            self.exitcode = EXIT_CODE_CONDA_ERROR if "Conda" in self.output else self.exitcode
            CONSOLE.print("STARTUP ERRORS")
        #sys.exit(0)

    def workflow_setup(self):
        while self.exitcode < 0:
            if self.quiet < 2:
                with CONSOLE.status("[dim]Preparing workflow", spinner = "point", spinner_style="yellow"):
                    while self.output.startswith("Building DAG of jobs...") or self.output.startswith("Assuming"):
                        self.nextline()
                if "Nothing to be" in self.output:
                    CONSOLE.rule("[bold]All outputs already present", style = "green")
                    sys.exit(0)
            else:
                while self.output.startswith("Building DAG of jobs...") or self.output.startswith("Assuming"):
                    self.nextline()
            if self.process.poll() or self.iserror():
                self.exitcode = EXIT_CODE_SUCCESS if self.process.poll() == 0 else EXIT_CODE_GENERIC_ERROR
                return
            while not self.output.startswith("Job stats:") and self.exitcode < 0:
                # print dependency text only once
                if "Downloading and installing remote packages" in self.output or "Running post-deploy" in self.output:
                    self.deps = True
                    self.deploy_text += "[dim]Installing workflow software"
                    break
                if "Pulling singularity image" in self.output:
                    self.deps = True
                    self.deploy_text += "[dim]Building software container"
                    break
                if "Nothing to be" in self.output:
                    CONSOLE.rule("[bold]All workflow outputs already present", style = "green")
                    sys.exit(0)
                if "MissingInput" in self.output:
                    self.exitcode = EXIT_CODE_GENERIC_ERROR
                    break
                if "AmbiguousRuleException" in self.output or "Error" in self.output or "Exception" in self.output:
                    self.exitcode = EXIT_CODE_RUNTIME_ERROR
                    break
                self.nextline()
            # if dependency text present, print pulsing progress bar
            if self.deps:
                progress = harpy_pulsebar(self.quiet)
                with harpy_progresspanel(progress, quiet=self.quiet, title = self.deploy_text):
                    progress.add_task("[dim]Working...", total = None)
                    while not self.output.startswith("Job stats:"):
                        self.nextline()
                        if self.process.poll() or self.iserror():
                            self.exitcode = EXIT_CODE_SUCCESS if self.process.poll() == 0 else 2
                            break
                    progress.stop()
            if self.process.poll() or self.exitcode >= 0:
                return
            if "Nothing to be" in self.output:
                CONSOLE.rule("[bold]All outputs already present", style = "green")
                sys.exit(0)
            # process the job summary
            while True:
                self.nextline()
                # stop parsing on "total" since it's always the last entry
                if self.output.startswith("Select jobs to execute"):
                    return
                try:
                    rule,count = self.output.split()
                    if rule in ["job", "all"] or "----" in rule:
                        continue
                    rule_desc = rule.replace("_", " ")
                    # rule : display_name, count_total, set of job_id's
                    self.job_inventory[rule] = [rule_desc, int(count), set()]
                except ValueError:
                    pass
            # checkpoint
                if self.process.poll() or self.iserror():
                    self.exitcode = EXIT_CODE_SUCCESS if self.process.poll() == 0 else EXIT_CODE_GENERIC_ERROR
                    return

    def monitor_jobs(self):
        self.job_inventory["total"][1] -= 1
        if self.exitcode > -1 or self.process.poll():
            return
        #CONSOLE.print("WE START MONITORING JOBS")
        progress = harpy_progressbar(self.quiet) 
        with harpy_progresspanel(progress, quiet = self.quiet):
            task_ids = {
                "total_progress" : progress.add_task(
                    "[bold blue]Total" if self.quiet == 0 else "[bold blue]Progress",
                    total= self.job_inventory["total"][1], active = 0
                )
            }

            while self.output:
                self.nextline()
                if self.iserror() or self.process.poll() == 1:
                    progress.stop()
                    #CONSOLE.log("WATATA")
                    self.exitcode = EXIT_CODE_RUNTIME_ERROR
                    break
                if "(100%) done" in self.output or self.output.startswith("Nothing to be") or self.process.poll() == 0:
                    progress.stop()
                    #CONSOLE.log("Alkaline")
                    self.exitcode = EXIT_CODE_SUCCESS
                    break
                if self.output.startswith("Complete log") or self.process.poll():
                    progress.stop()
                    #CONSOLE.log("RAINBOW")                    
                    self.exitcode = EXIT_CODE_SUCCESS if self.process.poll() == 0 else EXIT_CODE_RUNTIME_ERROR
                    break
                # add new progress bar track if the rule doesn't have one yet
                if self.output.lstrip().startswith("rule ") or self.output.lstrip().startswith("localrule "):
                    # catch the 2nd word and remove the colon
                    rule = self.output.split()[-1].replace(":", "")
                    # add progress bar if it doesn't exist
                    if rule not in task_ids and rule != "all":
                        task_ids[rule] = progress.add_task(self.job_inventory[rule][0], total=self.job_inventory[rule][1], visible = self.quiet != 1, active = 1)
                    # parse the rest of the rule block to get the job ID and add it to the inventory
                    while True:
                        self.nextline()
                        # store the job id in the inventory so we can later look up which rule it's associated with
                        if "jobid: " in self.output:
                            job_id = int(self.output.strip().split()[-1])
                            if rule != "all":
                                self.job_inventory[rule][2].add(job_id)
                                # update the number of active jobs
                                progress.update(task_ids[rule], active = len(self.job_inventory[rule][2]))
                                # update total
                                self.update_total_active()
                                progress.update(task_ids["total_progress"], refresh = True, active = self.total_active - 1)
                            break
                # check which rule the job is associated with and update the corresponding progress bar
                if self.output.startswith("Finished jobid: "):
                    completed = int(re.search(r"\d+", self.output).group())
                    for job,details in self.job_inventory.items():
                        if completed in details[2]:
                            task_id = task_ids[job]
                            current_task = progress.tasks[task_id]
                            progress.update(task_id, advance = 1, refresh = True, active = current_task.fields["active"] - 1)
                            self.update_total_active()
                            progress.update(task_ids["total_progress"], refresh=True, advance=1, active = self.total_active - 1)
                            if progress.tasks[task_ids[job]].completed == progress.tasks[task_id].total:
                                progress.update(task_ids[job], refresh=True, description=f"[dim]{details[0]}", active = "")
                                #progress.update(task_ids[job], visible = False)
                            # remove the job to save memory. wont be seen again
                            details[2].discard(completed)
                            break

    def process_finish(self):
        # keep going
        if self.exitcode == -1:
            return

        self.process.wait()
        if self.exitcode == 0:
            #print("PROCESS FINISH 0")
            return
        if self.exitcode in (1,2):
            print_setup_error(self.exitcode)
        elif self.exitcode == 3:
            print_onerror(os.path.join(os.path.basename(self.outdir), self.logfile), datetime.now() - self.start_time)
        CONSOLE.tab_size = 4
        CONSOLE._highlight = False

        while self.output:
            self.nextline()
            if "Exiting because a job execution failed" in self.output:
                self.nextline()
                if self.output.strip().startswith("[") and self.output.strip().endswith("]"):
                    # this is the [timestamp] line
                    CONSOLE.print("[blue]" + self.output.strip(), overflow = "ignore", crop = False)
                    break
        # error in rule line
        while self.output:
            self.nextline()
            if "Error in rule" in self.output:
                CONSOLE.print("[yellow bold]" + self.output.strip(), overflow = "ignore", crop = False)
            elif self.output.strip().startswith("shell:"):
                self.print_shellcmd()
            elif self.output.startswith("Complete log"):
                break
            elif self.output.startswith("WorkflowError"):
                break
            else:
                self.process_error()