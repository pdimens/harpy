import sys
import os
import re
import shutil
import subprocess
from .printfunctions import print_error, print_solution, print_solution_with_culprits
from collections import Counter
import rich_click as click

def createregions(infile, window, method):
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
        bedregion = []
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
                bedregion.append(f"{contig}:{startpos}-{endpos}")
        return bedregion, f"Genome/{bn}"

def generate_conda_deps():
    """Create the YAML files of the workflow conda dependencies"""
    condachannels = ["conda-forge", "bioconda", "defaults"]
    environ = {
        "qc" : ["falco", "fastp", "multiqc", "pysam=0.22"],
        "align": ["bwa", "ema","icu","libzlib", "minimap2", "samtools=1.19", "seqtk", "xz"],
        "variants.snp": ["bcftools=1.19", "freebayes=1.3.6"],
        "variants.sv": ["leviathan", "naibr-plus"],
        "phase" : ["hapcut2", "whatshap"],
        "simulations" : ["perl", "numpy"],
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

def fetch_file(file, destination, rename=None):
    """Find the 'file' in the PATH and copy it to the 'destination' with the original metadata"""
    result = None
    try:
        result = subprocess.check_output(["whereis", file])
    except:
        print_error(f"The GNU program \'whereis\', which is used to locate the file, was not found on the system and therefore unable to retrieve it. Terminating harpy.")
        print_solution("Make sure \'whereis\' is installed on the system. It is usually provided by default in all Unix-like operating systems.")

    if result is None:
        return []

    result = result.decode().splitlines()
    for line in result:
        if line.endswith(":"):
            print_error(f"The file \"{file}\" was not found in PATH, cannot run Harpy module.")
            print_solution(f"Make sure harpy was installed correctly and that you are in the harpy conda environment.")
            click.echo("See documentation: https://pdimens.github.io/harpy/install/")
            exit(1)
        else:
            result = line.split(" ")[-1]
            break
    
    os.makedirs(destination, exist_ok = True)
    if rename:
        destination += rename
    # copy2 to keep metadata during copy
    shutil.copy2(result, destination)

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
