"""launch snakemake"""

import re
import os
import sys
import glob
import subprocess
from datetime import datetime
from rich import print as rprint
from rich.console import Console
from ._misc import gzip_file, harpy_progressbar, harpy_pulsebar
from ._printing import print_onsuccess, print_onstart, print_onerror, print_setup_error

EXIT_CODE_SUCCESS = 0
EXIT_CODE_GENERIC_ERROR = 1
EXIT_CODE_CONDA_ERROR = 2
EXIT_CODE_RUNTIME_ERROR = 3
# quiet = 0 : print all things, full progressbar
# quiet = 1 : print all text, only "Total" progressbar
# quiet = 2 : print nothing, no progressbar
def iserror(text):
    """logical check for erroring trigger words in snakemake output"""
    return "Exception" in text or "Error" in text or "MissingOutputException" in text

def highlight_params(text):
    """make important snakemake attributes like 'input:' highlighted in the error output"""
    if text.startswith("    jobid:"):
        return text.replace("jobid:", "[bold default]jobid:[/]").rstrip()
    if text.startswith("    input:"):
        return text.replace("input:", "[bold default]input:[/]").rstrip()
    if text.startswith("    output:"):
        return text.replace("output:", "[bold default]output:[/]").rstrip()
    if text.startswith("    log:"):
        return text.replace("log:", "[bold default]log:[/]").rstrip()
    if text.startswith("    conda-env:"):
        return text.replace("conda-env:", "[bold default]conda-env:[/]").rstrip()
    if text.startswith("    container:"):
        return text.replace("container:", "[bold default]container:[/]").rstrip()
    if text.startswith("    shell:"): 
        return text.replace("shell:", "[bold default]shell:[/]").rstrip()
    if text.startswith("    wildcards:"): 
        return text.replace("wildcards:", "[bold default]wildcards:[/]").rstrip()
    if text.startswith("    affected files:"): 
        return text.replace("affected files:", "[bold default]affected files:[/]").rstrip()
    return text.rstrip()

def purge_empty_logs(target_dir):
    """scan target_dir and remove empty files, then scan it again and remove empty directories"""
    for logfile in glob.glob(f"{target_dir}/logs/**/*", recursive = True):
        if os.path.isfile(logfile) and os.path.getsize(logfile) == 0:
            os.remove(logfile)
    for logfile in glob.glob(f"{target_dir}/logs/**/*", recursive = True):
        if os.path.isdir(logfile) and not os.listdir(logfile):
            os.rmdir(logfile)

def launch_snakemake(sm_args, workflow, starttext, outdir, sm_logfile, quiet, summaryfile = None):
    """launch snakemake with the given commands"""
    start_time = datetime.now()
    if quiet < 2:
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
                exitcode = EXIT_CODE_SUCCESS if process.poll() == 0 else EXIT_CODE_GENERIC_ERROR
                exitcode = EXIT_CODE_CONDA_ERROR if "Conda" in output else exitcode
                break
            if quiet < 2:
                console = Console()
                with console.status("[dim]Preparing workflow", spinner = "point", spinner_style="yellow") as status:
                    while output.startswith("Building DAG of jobs...") or output.startswith("Assuming"):
                        output = process.stderr.readline()
                if "Nothing to be" in output:
                    exitcode = EXIT_CODE_SUCCESS
                    break
            else:
                while output.startswith("Building DAG of jobs...") or output.startswith("Assuming"):
                    output = process.stderr.readline()
            if process.poll() or iserror(output):
                exitcode = EXIT_CODE_SUCCESS if process.poll() == 0 else EXIT_CODE_GENERIC_ERROR
                break
            while not output.startswith("Job stats:") and not exitcode:
                # print dependency text only once
                if "Downloading and installing remote packages" in output or "Running post-deploy" in output:
                    deps = True
                    deploy_text = "[dim]Installing workflow software"
                    break
                if "Pulling singularity image" in output:
                    deps = True
                    deploy_text = "[dim]Building software container"
                    break
                if "Nothing to be" in output:
                    exitcode = EXIT_CODE_SUCCESS
                    break
                if "MissingInput" in output:
                    exitcode = EXIT_CODE_GENERIC_ERROR
                    break
                if "AmbiguousRuleException" in output or "Error" in output or "Exception" in output:
                    exitcode = EXIT_CODE_RUNTIME_ERROR
                    break
                output = process.stderr.readline()
            # if dependency text present, print pulsing progress bar
            if deps:
                with harpy_pulsebar(quiet, deploy_text) as progress:
                    progress.add_task("[dim]" + deploy_text, total = None)
                    while not output.startswith("Job stats:"):
                        output = process.stderr.readline()
                        if process.poll() or iserror(output):
                            exitcode = EXIT_CODE_SUCCESS if process.poll() == 0 else 2
                            break
                    progress.stop()
            if process.poll() or exitcode:
                break
            if "Nothing to be" in output:
                exitcode = EXIT_CODE_SUCCESS
                break
            with harpy_progressbar(quiet) as progress:
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
                        # rule : display_name, count_total, set of job_id's
                        job_inventory[rule] = [rule_desc, int(count), set()]
                    except ValueError:
                        pass
                # checkpoint
                    if process.poll() or iserror(output):
                        exitcode = EXIT_CODE_SUCCESS if process.poll() == 0 else EXIT_CODE_GENERIC_ERROR
                        break
                if process.poll() or exitcode:
                    break
                total_text = "[bold blue]Total" if quiet == 0 else "[bold blue]Progress"
                task_ids = {"total_progress" : progress.add_task(total_text, total=job_inventory["total"][1])}

                while output:
                    output = process.stderr.readline()
                    if iserror(output) or process.poll() == 1:
                        progress.stop()
                        exitcode = EXIT_CODE_RUNTIME_ERROR
                        break
                    if process.poll() == 0 or output.startswith("Complete log") or output.startswith("Nothing to be"):
                        progress.stop()
                        exitcode = EXIT_CODE_SUCCESS if process.poll() == 0 else EXIT_CODE_RUNTIME_ERROR
                        break
                    # add new progress bar track if the rule doesn't have one yet
                    #rulematch = re.search(r"(rule|checkpoint)\s\w+:", output)
                    if output.startswith("rule ") or output.startswith("localrule "):
                        # catch the 2nd word and remove the colon
                        rule = output.split()[-1].replace(":", "")
                        # add progress bar if it doesn't exist
                        if rule not in task_ids:
                            task_ids[rule] = progress.add_task(job_inventory[rule][0], total=job_inventory[rule][1], visible = quiet != 1)
                        # parse the rest of the rule block to get the job ID and add it to the inventory
                        while True:
                            output = process.stderr.readline()
                            if "jobid: " in output:
                                job_id = int(output.strip().split()[-1])
                                job_inventory[rule][2].add(job_id)
                                break
                        # store the job id in the inventory so we can later look up which rule it's associated with
                    # check which rule the job is associated with and update the corresponding progress bar
                    #finishmatch = re.search(r"Finished\sjobid:\s\d+", output)
                    if output.startswith("Finished jobid: "):
                        completed = int(re.search(r"\d+", output).group())
                        for job,details in job_inventory.items():
                            if completed in details[2]:
                                progress.advance(task_ids[job])
                                progress.advance(task_ids["total_progress"])
                                # remove the job to save memory. wont be seen again
                                details[2].discard(completed)
                                break

        process.wait()
        if process.returncode < 1:
            gzip_file(sm_logfile)
            purge_empty_logs(outdir)
            if quiet < 2:
                end_time = datetime.now()
                elapsed_time = end_time - start_time
                print_onsuccess(outdir, summaryfile, elapsed_time)
            sys.exit(0)
        else:
            if exitcode in (1,2):
                print_setup_error(exitcode)
            elif exitcode == 3:
                print_onerror(sm_logfile)
            while output and not output.endswith("]") and not output.startswith("Shutting down"):                   
                console = Console(stderr = True, tab_size = 4, highlight=False)
                if "Exception" in output or "Error" in output:
                    console.print("[yellow bold]" + output.rstrip(), overflow = "ignore", crop = False)
                    output = process.stderr.readline()
                    continue
                if output.startswith("Logfile") or output.startswith("======"):
                    console.print("[yellow bold]" + output.rstrip(), overflow = "ignore", crop = False)
                    output = process.stderr.readline()
                    continue
                if output:
                    if not output.startswith("Complete log"):
                        console.print("[red]" + highlight_params(output), overflow = "ignore", crop = False)
                    if output.startswith("Removing output files of failed job"):
                        break
                output = process.stderr.readline()
            gzip_file(sm_logfile)
            purge_empty_logs(outdir)
            sys.exit(1)
    except KeyboardInterrupt:
        # Handle a keyboard interrupt
        testfile.close()
        console = Console(stderr=True)
        console.print("")
        console.rule("[bold]Terminating Harpy", style = "yellow")
        process.terminate()
        process.wait()
        gzip_file(sm_logfile)
        purge_empty_logs(outdir)
        sys.exit(1)
