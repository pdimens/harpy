"""launch snakemake"""

import re
import sys
import subprocess
from rich import print as rprint
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, SpinnerColumn
from rich.console import Console
from ._printing import print_onsuccess, print_onstart, print_onerror, print_snakefile_error

def launch_snakemake(sm_args, workflow, starttext, outdir, sm_logfile, quiet):
    """launch snakemake with the given commands"""
    if not quiet:
        print_onstart(starttext, workflow.replace("_", " "))
    try:
        # Add a task with a total value of 100 (representing 100%)
        # Start a subprocess
        process = subprocess.Popen(sm_args.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
        err = False
        deps = False
        # read up to the job summary, but break early if dependency text appears
        while True:
            output = process.stderr.readline()
            return_code = process.poll()
            if return_code == 1:
                print_snakefile_error("If you manually edited the Snakefile, see the error below to fix it. If you didn't edit it manually, it's probably a bug (oops!) and you should submit an issue on GitHub: [bold]https://github.com/pdimens/harpy/issues")
                errtext = subprocess.run(sm_args.split(), stderr=subprocess.PIPE, text = True)
                errtext = errtext.stderr.split("\n")
                startprint = [i for i,j in enumerate(errtext) if "Exception in rule" in j or "Error in file" in j or "Exception:" in j or "Error:" in j or "Exception:" in j][0]
                for i in errtext[startprint:]:
                    if i:
                        rprint("[red]" + i)
                #rprint("[red]" + "\n".join(errtext[startprint:]), end = None)
                sys.exit(1)
            if not quiet:
                console = Console()
                with console.status("[dim]Preparing workflow", spinner = "point", spinner_style="yellow") as status:
                    while True:
                        if output.startswith("Building DAG of jobs...") or output.startswith("Assuming"):
                            pass
                        else:
                            break
                        output = process.stderr.readline()
            # print dependency text only once
            if "Downloading and installing remote packages" in output:
                deps = True
                deploy_text = "[dim]Installing software dependencies"
                break
            if "Pulling singularity image" in output:
                deps = True
                deploy_text = "[dim]Downloading software container"
                break
            if output.startswith("Job stats:"):
                break
        
        # if dependency text present, print pulse progress bar to indicate things are happening
        if deps:
            with Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(bar_width= 70 - len(deploy_text), pulse_style = "grey46"),
                TimeElapsedColumn(),
                transient=True,
                disable=quiet
            ) as progress:
                progress.add_task("[dim]" + deploy_text, total = None)
                while True:
                    output = process.stderr.readline()
                    if not output and process.poll():
                        progress.stop()
                        break
                    if output.startswith("Job stats:"):
                        progress.stop()
                        break
        with Progress(
            SpinnerColumn(spinner_name = "arc", style = "dim"),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(complete_style="yellow", finished_style="blue"),
            TextColumn("[progress.remaining]{task.completed}/{task.total}"),
            TimeElapsedColumn(),
            transient=True,
            disable=quiet
        ) as progress:
            # process the job summary
            job_inventory = {}
            while True:
                output = process.stderr.readline()
                if not output and process.poll():
                    break
                # stop parsing on "total" since it's always the last entry
                if output.startswith("total") or output.startswith("Select jobs to execute"):
                    break
                try:
                    rule,count = output.split()
                    rule_desc = rule.replace("_", " ")
                    # rule : display_name, count, set of job_id's
                    job_inventory[rule] = [rule_desc, int(count), set()]
                except ValueError:
                    pass
            task_ids = {"total_progress" : progress.add_task("[bold blue]Total", total=sum(j[1] for i,j in job_inventory.items()))}
            
            while True:
                output = process.stderr.readline()
                if process.poll():
                    if process.poll() == 1:
                        progress.stop()
                        print_snakefile_error("If you manually edited the Snakefile, see the error below to fix it. If you didn't edit it manually, it's probably a bug (oops!) and you should submit an issue on GitHub: [bold]https://github.com/pdimens/harpy/issues")
                        errtext = subprocess.run(sm_args.split(), stderr=subprocess.PIPE, text = True)
                        errtext = errtext.partition("total")
                        errtext = errtext.stderr.split("\n")
                        # if the DAG was created and job summary printed, print the stuff after
                        for i,j in enumerate(errtext):
                            if i.startswith("total"):
                                startfrom = j+1
                                for errline in errtext[startfrom:]:
                                    rprint("[red]" + errline)
                                sys.exit(1)
                        # otherwise, print everything
                        _ = [rprint("[red]" + i) for i in errtext]
                        sys.exit(1)
                    if not output:
                        break
                if output:
                    if output.startswith("Complete log:") or process.poll():
                        process.wait()
                        break
                    # prioritizing printing the error and quitting
                    if err:
                        while True:
                            output = process.stderr.readline()
                            if not output:
                                sys.exit(1)
                            if "Shutting down, this" in output or "%) done" in output or not output:
                                sys.exit(1)
                            if output.startswith("RuleException"):
                                rprint("[yellow bold]" + output.strip().partition("Finished job")[0], file = sys.stderr)
                            elif output.startswith("Error in rule"):
                                rprint("[yellow]" + output.strip().partition("Finished job")[0], file = sys.stderr)
                            else:
                                rprint("[red]" + output.strip().partition("Finished job")[0], file = sys.stderr)
                    if "Error in rule" in output or "RuleException" in output:
                        progress.stop()
                        if not err:
                            print_onerror(sm_logfile)
                        rprint(f"[yellow bold]{output.strip()}", file = sys.stderr)
                        err = True
                        continue
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
                                target_job = job
                                progress.update(task_ids[job], advance = 1)
                                progress.update(task_ids["total_progress"], advance=1)
                                # remove the job to save memory. wont be seen again
                                details[2].discard(completed)
                                break

        if process.returncode < 1:
            if not quiet:
                print_onsuccess(outdir)
        else:
            print_onerror(sm_logfile)
            with open(sm_logfile, "r", encoding="utf-8") as logfile:
                line = logfile.readline()
                while line:
                    if "Error" in line or "Exception" in line:
                        rprint("[bold yellow]" + line.rstrip())
                        break
                    line = logfile.readline()
                line = logfile.readline()
                while line:
                    if line.endswith("]\n"):
                        break
                    rprint("[red]" + line.rstrip())
                    line = logfile.readline()
            sys.exit(1)
    except KeyboardInterrupt:
        # Handle the keyboard interrupt
        rprint("[yellow bold]\nTerminating harpy...")
        process.terminate()
        process.wait()
        sys.exit(1)
