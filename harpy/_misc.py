"""Module with helper function to set up Harpy workflows"""

import os
import sys
import glob
import click
import subprocess
from datetime import datetime
from pathlib import Path
from collections import Counter
from importlib_resources import files
import harpy.scripts
import harpy.reports
import harpy.snakefiles
from ._printing import print_error, print_solution

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


class IntPair(click.ParamType):
    """A class for a click type which accepts 2 integers, separated by a comma."""
    name = "int_pair"
    def convert(self, value, param, ctx):
        try:
            parts = value.split(',')
            if len(parts) != 2:
                self.fail(f"{value} is not a valid int pair. The value should be two integers separated by a comma.", param, ctx)
            return [int(i) for i in parts]
        except ValueError:
            self.fail(f"{value} is not a valid int pair. The value should be two integers separated by a comma.", param, ctx)

class IntQuartet(click.ParamType):
    """A class for a click type which accepts 4 integers, separated by a comma."""
    name = "int_pair"
    def convert(self, value, param, ctx):
        try:
            parts = value.split(',')
            if len(parts) != 4:
                self.fail(f"{value} is not a valid set of 4 integers separated by a comma.", param, ctx)
            return [int(i) for i in parts]
        except ValueError:
            self.fail(f"{value} is not a valid set of 4 integers separated by a comma.", param, ctx)
