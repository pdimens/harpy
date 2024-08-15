"""Module with helper function to set up Harpy workflows"""

import os
import re
import sys
import glob
import subprocess
from time import sleep
from datetime import datetime
from pathlib import Path
from collections import Counter
from rich import print as rprint
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn
from rich.console import Console
from importlib_resources import files
import harpy.scripts
import harpy.reports
import harpy.snakefiles
from .printfunctions import print_error, print_solution, print_onsuccess, print_onstart, print_onerror, print_snakefile_error

def symlink(original, destination):
    """Create a symbolic link from original -> destination if the destination doesn't already exist."""
    if not (Path(destination).is_symlink() or Path(destination).exists()):
        Path(destination).symlink_to(Path(original).resolve())

def fetch_script(workdir, target):
    """
    Retrieve the target harpy script and write it into workdir/scripts
    """
    os.makedirs(f"{workdir}/scripts/", exist_ok= True)
    with open(f"{workdir}/scripts/{target}", "w", encoding="utf-8") as f:
        if os.path.isfile(files(harpy.scripts).joinpath(target)):
            f.write(files(harpy.scripts).joinpath(target).read_text())
        else:
            print_error("script missing", f"Bundled script [blue bold]{target}[/blue bold] was not found in the Harpy installation.")
            print_solution("There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration.")
            sys.exit(1)

def fetch_rule(workdir, target):
    """
    Retrieve the target harpy rule and write it into the workdir
    """
    os.makedirs(f"{workdir}/", exist_ok= True)
    with open(f"{workdir}/{target}", "w", encoding="utf-8") as f:
        if os.path.isfile(files(harpy.snakefiles).joinpath(target)):
            f.write(files(harpy.snakefiles).joinpath(target).read_text())
        else:
            print_error("snakefile missing", f"The required snakefile [blue bold]{target}[/blue bold] was not found in the Harpy installation.")
            print_solution("There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration.")
            sys.exit(1)

def fetch_report(workdir, target):
    """
    Retrieve the target harpy report and write it into workdir/report
    """
    os.makedirs(f"{workdir}/report/", exist_ok= True)
    with open(f"{workdir}/report/{target}", "w", encoding="utf-8") as f:
        if os.path.isfile(files(harpy.reports).joinpath(target)):
            f.write(files(harpy.reports).joinpath(target).read_text())
        else:
            print_error("report script missing", f"The required report script [blue bold]{target}[/blue bold] was not found within the Harpy installation.")
            print_solution("There may be an issue with your Harpy installation, which would require reinstalling Harpy. Alternatively, there may be in a issue with your conda/mamba environment or configuration.")
            sys.exit(1)

def biallelic_contigs(vcf, workdir):
    """Identify which contigs have at least 2 biallelic SNPs"""
    vbn = os.path.basename(vcf)
    if not os.path.exists(f"{workdir}/{vbn}.biallelic"):
        os.makedirs(f"{workdir}/", exist_ok = True)
        biallelic = subprocess.Popen(f"bcftools view -M2 -v snps {vcf} -Ob".split(), stdout = subprocess.PIPE)
        contigs = subprocess.run("""bcftools query -f '%CHROM\\n'""".split(), stdin = biallelic.stdout, stdout = subprocess.PIPE, check = False).stdout.decode().splitlines()
        counts = Counter(contigs)
        contigs = [i.replace("\'", "") for i in counts if counts[i] > 1]
        with open(f"{workdir}/{vbn}.biallelic", "w", encoding="utf-8") as f:
            _ = [f.write(f"{i}\n") for i in contigs]
    else:
        with open(f"{workdir}/{vbn}.biallelic", "r", encoding="utf-8") as f:
            contigs = [line.rstrip() for line in f]
    if len(contigs) == 0:
        print_error("no usable contigs", "No contigs with at least 2 biallelic SNPs identified. Cannot continue with imputation.")
        sys.exit(1)
    else:
        return contigs

def snakemake_log(outdir, workflow):
    """Return a snakemake logfile name. Iterates logfile run number if one exists."""
    attempts = glob.glob(f"{outdir}/logs/snakemake/*.log")
    if not attempts:
        return f"{outdir}/logs/snakemake/{workflow}.1." + datetime.now().strftime("%d_%m_%Y") + ".log"
    increment = sorted([int(i.split(".")[1]) for i in attempts])[-1] + 1
    return f"{outdir}/logs/snakemake/{workflow}.{increment}." + datetime.now().strftime("%d_%m_%Y") + ".log"

def launch_snakemake(sm_args, workflow, starttext, outdir, sm_logfile, quiet):
    """launch snakemake with the given commands"""
    print_onstart(starttext, workflow.replace("_", " "))
    try:
        with Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            transient=True,
            disable=quiet
        ) as progress:
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
                    rprint("[red]" + errtext.stderr.partition("jobs...\n")[2], end = None)
                    sys.exit(1)
                if output == '' and return_code is not None:
                    break
                # print dependency text only once
                if "Downloading and installing remote packages" in output:
                    deps = True
                    deploy_text = "[magenta]Downloading and installing workflow dependencies"
                    break
                if "Pulling singularity image" in output:
                    deps = True
                    deploy_text = "[magenta]Downloading software container"
                    break
                if output.startswith("Job stats:"):
                    # read and ignore the next two lines
                    process.stderr.readline()
                    process.stderr.readline()
                    break
            # if dependency text present, print console log with spinner and read up to the job stats
            if deps:
                if not quiet:
                    console = Console()
                    with console.status(deploy_text, spinner = "point") as status:
                        while True:
                            output = process.stderr.readline()
                            if output == '' and process.poll() is not None:
                                break
                            sleep(2)
                            if output.startswith("Job stats:"):
                                # read and ignore the next two lines
                                process.stderr.readline()
                                process.stderr.readline()
                                break
                else:
                    while True:
                        output = process.stderr.readline()
                        if output == '' and process.poll() is not None:
                            break
                        if output.startswith("Job stats:"):
                            # read and ignore the next two lines
                            process.stderr.readline()
                            process.stderr.readline()
                            break
            job_inventory = {}
            task_ids = {"total_progress" : progress.add_task("[bold blue]Total", total=100)}
            # process the job summary
            while True:
                output = process.stderr.readline()
                # stop parsing on "total" since it's always the last entry
                if output.startswith("total"):
                    break
                if output.startswith("workflow_summary"):
                    continue
                try:
                    rule,count = output.split()
                    rule_desc = rule.replace("_", " ")
                    # rule : display_name, count, set of job_id's
                    job_inventory[rule] = [rule_desc, int(count), set()]
                except ValueError:
                    pass

            while True:
                output = process.stderr.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    # prioritizing printing the error and quitting
                    if err:
                        if "Shutting down, this" in output:
                            sys.exit(1)
                        rprint(f"[red]{output.strip()}", file = sys.stderr)
                    if "Error in rule" in output:
                        progress.stop()
                        print_onerror(sm_logfile)
                        rprint(f"[yellow bold]{output.strip()}", file = sys.stderr)
                        err = True
                    # find the % progress text to update Total progress bar
                    match = re.search(r"\(\d+%\) done", output)
                    if match:
                        percent = int(re.sub(r'\D', '', match.group()))
                        progress.update(task_ids["total_progress"], completed=percent)
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
                                progress.update(task_ids[job], advance = 1)
                                # remove the job to save memory. wont be seen again
                                details[2].discard(completed)
                                break

        if process.returncode < 1:
            print_onsuccess(outdir)
    except KeyboardInterrupt:
        # Handle the keyboard interrupt
        rprint("[yellow bold]\nTerminating harpy...")
        process.terminate()
        process.wait()
        sys.exit(1)
