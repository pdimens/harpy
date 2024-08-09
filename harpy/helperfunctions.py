"""Module with helper function to set up Harpy workflows"""

import os
import re
import sys
import subprocess
from pathlib import Path
from collections import Counter
from rich import print as rprint
from rich.progress import Progress
from importlib_resources import files
import harpy.scripts
import harpy.reports
import harpy.snakefiles
from .printfunctions import print_error, print_solution, print_success

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
            print_error(f"Bundled script [blue bold]{target}[/blue bold] was not found in the Harpy installation.")
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
            print_error(f"The required snakefile [blue bold]{target}[/blue bold] was not found in the Harpy installation.")
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
            print_error(f"The required report script [blue bold]{target}[/blue bold] was not found within the Harpy installation.")
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
        print_error("No contigs with at least 2 biallelic SNPs identified. Cannot continue with imputation.")
        sys.exit(1)
    else:
        return contigs

def launch_snakemake(sm_args, outdir):
    """launch snakemake with the given commands"""
    try:
        with Progress(transient=True) as progress:
        # Add a task with a total value of 100 (representing 100%)
            task = progress.add_task("harpy qc ", total=100)
            # Start a subprocess
            process = subprocess.Popen(sm_args.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
            err = 0
            while True:
                output = process.stderr.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    if err > 0:
                        rprint(f"[red]{output.strip()}")
                    if "Error in rule" in output:
                        #TODO PRINT ERROR BOX HERE
                        rprint(f"[yellow bold]{output.strip()}")
                        err += 1
                    match = re.search(r"\(\d+%\)", output)
                    if match:
                        percent = int(re.sub(r'\D', '', match.group()))
                        progress.update(task, completed=percent)
            print_success(outdir)
    except KeyboardInterrupt:
        # Handle the keyboard interrupt
        rprint("[yellow bold]\nTerminating harpy...")
        process.terminate()
        process.wait()
        sys.exit(1)

#TODO PATCH THIS UP, ADD IT TO FUNCTIONS
def new_snakemake_logfile():
    attempts = glob.glob(f"{outdir}/logs/snakemake/*.snakelog")
    if not attempts:
        logfile = f"{outdir}/logs/snakemake/qc.run1." + datetime.now().strftime("%d_%m_%Y") + ".snakelog"
    else:
        increment = sorted([int(i.split(".")[1].replace("run","")) for i in attempts])[-1] + 1
        logfile = f"{outdir}/logs/snakemake/qc.run{increment}." + datetime.now().strftime("%d_%m_%Y") + ".snakelog"
