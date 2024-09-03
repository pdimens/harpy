"""launch snakemake"""

import re
import sys
import subprocess
from rich import print as rprint
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, SpinnerColumn
from rich.console import Console
from ._misc import gzip_file
from ._printing import print_onsuccess, print_onstart, print_onerror, print_setup_error

def iserror(text):
    """logical check for erroring trigger words in snakemake output"""
    return "Exception" in text or "Error" in text or "MissingOutputException" in text

def launch_snakemake(sm_args, workflow, starttext, outdir, sm_logfile, quiet):
    """launch snakemake with the given commands"""
    if not quiet:
        print_onstart(starttext, workflow.replace("_", " "))
    exitcode = None
    try:
        # Start snakemake as a subprocess
        process = subprocess.Popen(sm_args.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
        deps = False
        # read up to the job summary, but break early if dependency text appears
        while not exitcode:
            output = process.stderr.readline()
            # check for syntax errors at the very beginning
            if process.poll() or iserror(output):
                exitcode = 0 if process.poll() == 0 else 1
                break
            if not quiet:
                console = Console()
                with console.status("[dim]Preparing workflow", spinner = "point", spinner_style="yellow") as status:
                    while output.startswith("Building DAG of jobs...") or output.startswith("Assuming"):
                        output = process.stderr.readline()
            else:
                while output.startswith("Building DAG of jobs...") or output.startswith("Assuming"):
                    output = process.stderr.readline()
            if process.poll() or iserror(output):
                exitcode = 0 if process.poll() == 0 else 1

            while not output.startswith("Job stats:"):
                # print dependency text only once
                if "Downloading and installing remote packages" in output:
                    deps = True
                    deploy_text = "[dim]Installing software dependencies"
                    break
                if "Pulling singularity image" in output:
                    deps = True
                    deploy_text = "[dim]Downloading software container"
                    break
                output = process.stderr.readline()

            # if dependency text present, print pulsing progress bar
            if deps:
                with Progress(
                    TextColumn("[progress.description]{task.description}"),
                    BarColumn(bar_width= 70 - len(deploy_text), pulse_style = "grey46"),
                    TimeElapsedColumn(),
                    transient=True,
                    disable=quiet
                ) as progress:
                    progress.add_task("[dim]" + deploy_text, total = None)
                    while not output.startswith("Job stats:"):
                        output = process.stderr.readline()
                        if process.poll() or iserror(output):
                            exitcode = 0 if process.poll() == 0 else 2
                            break
                    progress.stop()
            if process.poll() or exitcode:
                break
            with Progress(
                SpinnerColumn(spinner_name = "arc", style = "dim"),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(complete_style="yellow", finished_style="blue"),
                TextColumn("[progress.remaining]{task.completed}/{task.total}", style = "magenta"),
                TimeElapsedColumn(),
                transient=True,
                disable=quiet
            ) as progress:
                # process the job summary
                job_inventory = {}
                while True:
                    output = process.stderr.readline()
                    # stop parsing on "total" since it's always the last entry
                    if output.startswith("Select jobs to execute"):
                        break
                    try:
                        rule,count = output.split()
                        rule_desc = rule.replace("_", " ")
                        # rule : display_name, count, set of job_id's
                        job_inventory[rule] = [rule_desc, int(count), set()]
                    except ValueError:
                        pass
                # checkpoint
                    if process.poll() or iserror(output):
                        exitcode = 0 if process.poll() == 0 else 1
                        break
                if process.poll() or exitcode:
                    break
                task_ids = {"total_progress" : progress.add_task("[bold blue]Total", total=job_inventory["total"][1])}

                while output:
                    output = process.stderr.readline()
                    if iserror(output) or process.poll() == 1:
                        progress.stop()
                        exitcode = 3
                        break
                    if process.poll() == 0 or output.startswith("Complete log:"):
                        progress.stop()
                        exitcode = 0 if process.poll() == 0 else 3
                        break
                    # add new progress bar track if the rule doesn't have one yet
                    rulematch = re.search(r"rule\s\w+:", output)
                    if rulematch:
                        rule = rulematch.group().replace(":","").split()[-1]
                        if rule not in task_ids:
                            task_ids[rule] = progress.add_task(job_inventory[rule][0], total=job_inventory[rule][1])
                        continue
                    # store the job id in the inventory so we can later look up which rule it's associated with
                    jobidmatch = re.search(r"jobid:\s\d+", string = output)
                    if jobidmatch:
                        job_id = int(re.search(r"\d+",output).group())
                        # rule should be the most previous rule recorded
                        job_inventory[rule][2].add(job_id)
                        continue
                    # check which rule the job is associated with and update the corresponding progress bar
                    finishmatch = re.search(r"Finished\sjob\s\d+", output)
                    if finishmatch:
                        completed = int(re.search(r"\d+", output).group())
                        for job,details in job_inventory.items():
                            if completed in details[2]:
                                progress.update(task_ids[job], advance = 1)
                                progress.update(task_ids["total_progress"], advance=1)
                                # remove the job to save memory. wont be seen again
                                details[2].discard(completed)
                                break

        process.wait()
        if process.returncode < 1:
            gzip_file(sm_logfile)
            if not quiet:
                print_onsuccess(outdir)
            sys.exit(0)
        else:
            if exitcode in (1,2):
                print_setup_error(exitcode)
            elif exitcode == 3:
                print_onerror(sm_logfile)
            while output and not output.endswith("]") and not output.startswith("Shutting down"):
                if "Exception" in output or "Error" in output:
                    rprint("[yellow bold]" + output.rstrip(), file = sys.stderr)
                    output = process.stderr.readline()
                    continue
                if output.startswith("Logfile"):
                    rprint("[yellow]" + output.rstrip(), file = sys.stderr)
                    output = process.stderr.readline()
                    continue
                if output:
                    if not output.startswith("Complete log"):
                        rprint("[red]" + output.replace("\t","    ").rstrip(), file = sys.stderr)
                output = process.stderr.readline()
            sys.exit(1)
    except KeyboardInterrupt:
        # Handle the keyboard interrupt
        rprint("[yellow bold]\nTerminating harpy...", file = sys.stderr)
        process.terminate()
        process.wait()
        gzip_file(sm_logfile)
        sys.exit(1)
