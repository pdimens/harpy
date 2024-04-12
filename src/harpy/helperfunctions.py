import sys
import os
import re
import shutil
import subprocess
from pathlib import Path
from importlib_resources import files
from .printfunctions import print_error, print_solution, print_solution_with_culprits
from collections import Counter
import rich_click as click
import harpy.scripts
import harpy.reports
import harpy.rules

def symlink(original, destination):
    """Create a symbolic link from original -> destination if the destination doesn't already exist."""
    if not (Path(destination).is_symlink() or Path(destination).exists()):
        Path(destination).symlink_to(Path(original).absolute()) 

def createregions(infile, window, method, outfile):
    """Create a BED file of genomic intervals of size 'window'. Uses 1- or 0- based numbering depending on mpileup or freebayes 'method'"""
    bn = os.path.basename(infile)
    os.makedirs("Genome", exist_ok = True)
    base = 0 if method == "freebayes" else 1
    gen_zip = True if bn.lower().endswith(".gz") else False
    if method == "freebayes":
        # freebayes requires uncompressed genome
        if gen_zip:
            # remove .gz extension
            bn = bn[:-3]
            if not os.path.exists(f"Genome/{bn}"):
                with open(f"Genome/{bn}", "w") as fo:
                    subprocess.run(f"gzip -dc {infile}".split(), stdout = fo)
        else:
            if not os.path.exists(f"Genome/{bn}"):
                subprocess.run(f"ln -sr {infile} Genome/{bn}".split())
    else:
        if not os.path.exists(f"Genome/{bn}"):
            ftype = subprocess.run(["file", infile], stdout=subprocess.PIPE).stdout.decode('utf-8')
            if "Blocked GNU Zip" in ftype:
                # is bgzipped, just link it
                subprocess.run(f"ln -sr {infile} Genome/{bn}".split())
            elif "gzip compressed data" in ftype:
                # is regular gzipped, needs to be bgzipped
                subprocess.run(f"zcat {infile} | bgzip -c > Genome/{bn}".split())
            else:
                # not compressed, just link
                subprocess.run(f"ln -sr {infile} Genome/{bn}".split())

    if not os.path.exists(f"Genome/{bn}.fai"):
        try:
            subprocess.run(f"samtools faidx --fai-idx Genome/{bn}.fai --gzi-idx Genome/{bn}.gzi Genome/{bn}".split(), stderr = subprocess.DEVNULL)
        except:
            subprocess.run(f"samtools faidx --fai-idx Genome/{bn}.fai Genome/{bn}".split(), stderr = subprocess.DEVNULL)

    with open(f"Genome/{bn}.fai") as fai:
        bedout = open(outfile, "w")
        while True:
            # Get next line from file
            line = fai.readline()
            # if line is empty, end of file is reached
            if not line:
                break
            # split the line by tabs
            lsplit = line.split()
            contig = lsplit[0]
            c_len = int(lsplit[1])
            c_len = c_len - 1 if base == 0 else c_len
            start = base
            end = window
            starts = [base]
            ends = [window]
            while end < c_len:
                end = end + window if (end + window) < c_len else c_len
                ends.append(end)
                start += window
                starts.append(start)
            for (startpos, endpos) in zip (starts,ends):
                bedout.write(f"{contig}:{startpos}-{endpos}\n")
        bedout.close()
        return f"Genome/{bn}"

def generate_conda_deps():
    """Create the YAML files of the workflow conda dependencies"""
    condachannels = ["conda-forge", "bioconda", "defaults"]
    environ = {
        "qc" : ["falco", "fastp", "multiqc", "pysam=0.22"],
        "align": ["bwa", "ema","icu","libzlib", "minimap2", "samtools=1.19", "seqtk", "xz"],
        "variants.snp": ["bcftools=1.19", "freebayes=1.3.6"],
        "variants.sv": ["leviathan", "naibr-plus"],
        "phase" : ["hapcut2", "whatshap"],
        "simulations" : ["perl", "perl-math-random", "perl-inline-c", "perl-parse-recdescent", "numpy", "dwgsim", "alienzj::msort"],
        "r-env" : ["bioconductor-complexheatmap", "r-circlize", "r-dt", "r-flexdashboard", "r-ggplot2", "r-ggridges", "r-plotly", "r-tidyr", "r-stitch"]
    }

    os.makedirs(".harpy_envs", exist_ok = True)

    for i in environ:
        # don't overwrite existing
        if not os.path.isfile(f".harpy_envs/{i}.yaml"):
            with open(f".harpy_envs/{i}.yaml", "w") as yml:
                yml.write(f"name: {i}\n")
                yml.write("channels:\n  - ")
                yml.write("\n  - ".join(condachannels))
                yml.write("\ndependencies:\n  - ")
                yml.write("\n  - ".join(environ[i]) + "\n")

def fetch_script(workdir, target):
    """
    Retrieve the target harpy script and write it into workdir/scripts
    """
    os.makedirs(f"{workdir}/scripts/", exist_ok= True)
    with open(f"{workdir}/scripts/{target}", "w") as f:
        f.write(files(harpy.scripts).joinpath(target).read_text())

def fetch_rule(workdir, target):
    """
    Retrieve the target harpy rule and write it into the workdir
    """
    os.makedirs(f"{workdir}/", exist_ok= True)
    with open(f"{workdir}/{target}", "w") as f:
        f.write(files(harpy.rules).joinpath(target).read_text())

def fetch_report(workdir, target):
    """
    Retrieve the target harpy report and write it into workdir/report
    """
    os.makedirs(f"{workdir}/report/", exist_ok= True)
    with open(f"{workdir}/report/{target}", "w") as f:
        f.write(files(harpy.reports).joinpath(target).read_text())

def biallelic_contigs(vcf):
    """Identify which contigs have at least 2 biallelic SNPs"""
    vbn = os.path.basename(vcf)
    if not os.path.exists(f"Impute/input/_{vbn}.list"):
        os.makedirs("Impute/input/", exist_ok = True)
        biallelic = subprocess.Popen(f"bcftools view -M2 -v snps {vcf} -Ob".split(), stdout = subprocess.PIPE)
        contigs = subprocess.run("""bcftools query -f '%CHROM\\n'""".split(), stdin = biallelic.stdout, stdout = subprocess.PIPE).stdout.decode().splitlines()
        counts = Counter(contigs)
        contigs = [i.replace("\'", "") for i in counts if counts[i] > 1]
        with open(f"Impute/input/_{vbn}.list", "w") as f:
            _ = [f.write(f"{i}\n") for i in contigs]
    else:
        with open(f"Impute/input/_{vbn}.list", "r") as f:
            contigs = [line.rstrip() for line in f]
    
    if len(contigs) == 0:
        print_error("No contigs with at least 2 biallelic SNPs identified. Cannot continue with imputation.")
        exit(1)
    else:
        return contigs
