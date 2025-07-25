"""Module for workflow file parsers"""

import re
import os
import sys
import subprocess
from typing import Tuple
from pathlib import Path
from rich.markdown import Markdown
import rich_click as click
from .misc import filepath
from .printing import print_error, print_solution_with_culprits

def getnames(directory: str, ext: str) -> list[str]:
    """Find all files in 'directory' that end with 'ext'"""
    samplenames = set([i.split(ext)[0] for i in os.listdir(directory) if i.endswith(ext)])
    if len(samplenames) < 1:
        print_error("no files found", f"No sample files ending with [bold]{ext}[/] found in [bold]{directory}[/].")
        sys.exit(1)
    return samplenames

def parse_impute_regions(regioninput: str, vcf: str) -> list:
    """Validates the --regions input of harpy impute to infer it's a properly formatted region
    Use the contigs and lengths of the vcf file to check that the region is valid. Returns
    a tuple of (contig, start, end)
    """
    contigs = contigs_from_vcf(vcf)
    # already validated by CLI type checking
    contig, positions = regioninput.split(":")
    startpos,endpos,buffer = [int(i) for i in positions.split("-")]
    # check if the region is in the genome
    if contig not in contigs:
        print_error("contig not found", f"The contig [bold yellow]{contig}[/] was not found in [blue]{vcf}[/].")
        sys.exit(1)
    if endpos > contigs[contig]:
        print_error("invalid region", f"The region end position [yellow bold]({endpos})[/] is greater than the length of contig [yellow bold]{contig}[/] ({contigs[contig]})")
        sys.exit(1)
    return contig, startpos, endpos, buffer

def parse_fastq_inputs(inputs: list[str], param: str = None) -> Tuple[list[str], int]:
    """
    Parse the command line input FASTQ arguments to generate a clean list of input files. Returns the number of unique samples,
    i.e. forward and reverse reads for one sample = 1 sample.
    """
    infiles = []
    re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
    for i in inputs:
        if os.path.isdir(i):
            for j in os.listdir(i):
                if re.search(re_ext, j):
                    infiles.append(filepath(os.path.join(i, j)))
        else:
            if re.search(re_ext, i):
                infiles.append(filepath(i))
    if len(infiles) < 1:
        print_error("no fastq files found", f"There were no files ending with the accepted fastq extensions [blue].fq[dim][.gz][/][/] or [blue].fastq[dim][.gz][/][/] in the provided [green]{param}[/] (case insensitive).")
        sys.exit(1)
    # check if any names will be clashing
    bn_r = r"[\.\_](?:[RF])?(?:[12])?(?:\_00[1-9])*?$"
    uniqs = set()
    dupes = []
    inv_pattern = r'[^a-zA-Z0-9._-]+'
    badmatch = []
    for i in infiles:
        sans_ext = os.path.basename(re_ext.sub("", str(i)))
        if sans_ext in uniqs:
            dupes.append(sans_ext)
        else:
            uniqs.add(sans_ext)
        if re.search(inv_pattern, os.path.basename(i)):
            badmatch.append(os.path.basename(i))
    if badmatch:
        print_error("invalid characters", f"Invalid characters were detected in the file names for [green]{param}[/].")
        print_solution_with_culprits(Markdown("Valid file names may contain only:\n- **A-Z** characters (case insensitive)\n- **.** (period)\n- **_** (underscore)\n- **-** (dash)"), "The offending files:")
        click.echo(", ".join(badmatch), file = sys.stderr)
        sys.exit(1)
    if dupes:
        print_error("clashing sample names", Markdown("Identical filenames were detected in `{param}`, which will cause unexpected behavior and results.\n- files with identical names but different-cased extensions are treated as identical\n- files with the same name from different directories are also considered identical"))
        print_solution_with_culprits("Make sure all input files have unique names.", "Files with clashing names:")
        for i in dupes:
            click.echo(" ".join([j for j in infiles if i in j]), file = sys.stderr)
        sys.exit(1)

    n = len({re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in uniqs})
    # return the filenames and # of unique samplenames
    return infiles, n

def parse_alignment_inputs(inputs:list[str], param: str = None) -> Tuple[list[str], int]:
    """
    Parse the command line input sam/bam arguments to generate a clean list of input files
    and return the number of unique samples.
    """
    bam_infiles = []
    bai_infiles = []
    re_bam = re.compile(r".*\.(bam|sam)$", flags = re.IGNORECASE)
    re_bai = re.compile(r".*\.bam\.bai$", flags = re.IGNORECASE)
    for i in inputs:
        if os.path.isdir(i):
            for j in os.listdir(i):
                if re_bam.match(j):
                    bam_infiles.append(filepath(os.path.join(i, j)))
                elif re_bai.match(j):
                    bai_infiles.append(filepath(os.path.join(i, j)))
        else:
            if re_bam.match(i):
                bam_infiles.append(filepath(i))
            elif re_bai.match(i):
                bai_infiles.append(filepath(i))
    if len(bam_infiles) < 1:
        print_error("no bam files found", "There were no files ending with the [blue].bam[/] or [blue].sam[/] extensions in the provided {param} (case insensitive).")
        sys.exit(1)
    re_ext = re.compile(r"\.(bam|sam)$", re.IGNORECASE)

    # check if any links will be clashing
    uniqs = set()
    dupes = []
    inv_pattern = r'[^a-zA-Z0-9._-]+'
    badmatch = []
    for i in bam_infiles:
        bn = os.path.basename(re_ext.sub("", str(i)))
        if bn in uniqs:
            dupes.append(bn)
        else:
            uniqs.add(bn)
        if re.search(inv_pattern, os.path.basename(i)):
            badmatch.append(os.path.basename(i))
    if badmatch:
        print_error("invalid characters", "Invalid characters were detected in the input file names.")
        print_solution_with_culprits(Markdown("Valid file names may contain only:\n- **A-Z** characters (case insensitive)\n- **.** (period)\n- **_** (underscore)\n- **-** (dash)"), "The offending files:")
        click.echo(", ".join(badmatch), file = sys.stderr)
        sys.exit(1)
    if dupes:
        print_error("clashing sample names", Markdown("Identical filenames were detected, which will cause unexpected behavior and results.\n- files with identical names but different-cased extensions are treated as identical\n- files with the same name from different directories are also considered identical"))
        print_solution_with_culprits("Make sure all input files have unique names.", "Files with clashing names:")
        for i in dupes:
            click.echo(" ".join([j for j in bam_infiles if i in j]), file = sys.stderr)
        sys.exit(1)
    return bam_infiles, len(uniqs)

def biallelic_contigs(vcf: str, workdir: str) -> Tuple[str,list[str], int]:
    """Identify which contigs have at least 2 biallelic SNPs and write them to workdir/vcf.biallelic"""
    vbn = os.path.basename(vcf)
    os.makedirs(workdir, exist_ok = True)
    valid = []
    vcfheader = subprocess.check_output(['bcftools', 'view', '-h', vcf]).decode().split('\n')
    header_contigs = [i.split(",")[0].replace("##contig=<ID=","") for i in vcfheader if i.startswith("##contig=")]
    if vcf.lower().endswith("bcf") and not os.path.exists(f"{vcf}.csi"):
        subprocess.run(f"bcftools index {vcf}".split())
    if vcf.lower().endswith("vcf.gz") and not os.path.exists(f"{vcf}.tbi"):
        subprocess.run(f"bcftools index --tbi {vcf}".split())

    for contig in header_contigs:
        # Use bcftools to count the number of biallelic SNPs in the contig
        viewcmd = subprocess.Popen(['bcftools', 'view', '-H', '-r', contig, '-v', 'snps', '-m2', '-M2', '-c', '2', vcf], stdout=subprocess.PIPE)
        snpcount = 0
        while True:
            # Read the next line of output
            line = viewcmd.stdout.readline().decode()
            if not line:
                break
            snpcount += 1
            # If there are at least 2 biallellic snps, terminate the process
            if snpcount >= 2:
                valid.append(contig)
                viewcmd.terminate()
                break
    if not valid:
        click.echo("No contigs with at least 2 biallelic SNPs identified. Cannot continue with imputation.")
        sys.exit(1)
    with open(f"{workdir}/{vbn}.biallelic", "w", encoding="utf-8") as f:
        f.write("\n".join(valid))
    return filepath(f"{workdir}/{vbn}.biallelic"), valid, len(valid)

def contigs_from_vcf(vcf: str) -> dict:
    """reads the header of a vcf/bcf file and returns a dict of the contigs (keys) and their lengths (values)"""
    header = subprocess.check_output(f"bcftools head {vcf}".split(), text = True)
    contigs = [i for i in header.splitlines() if i.startswith("##contig=<ID")]
    d = {}
    for i in contigs:
        contig,length = i.split(",length=")
        contig = contig.replace("##contig=<ID=", "")
        # remove the trailing '>' and convert to integer
        length = int(length[:-1])
        d[contig] = length
    return d