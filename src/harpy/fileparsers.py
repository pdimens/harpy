"""Module for workflow file parsers"""

import re
import os
import sys
import subprocess
from pathlib import Path
from rich.markdown import Markdown
import rich_click as click
from .printfunctions import print_error, print_solution_with_culprits

def getnames(directory, ext):
    """Find all files in 'directory' that end with 'ext'"""
    samplenames = set([i.split(ext)[0] for i in os.listdir(directory) if i.endswith(ext)])
    if len(samplenames) < 1:
        print_error(f"No sample files ending with [bold]{ext}[/bold] found in [bold]{directory}[/bold].")
        sys.exit(1)
    return samplenames

def parse_fastq_inputs(inputs):
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
                    infiles.append(Path(os.path.join(i, j)).resolve())
        else:
            if re.search(re_ext, i):
                infiles.append(Path(i).resolve())

    if len(infiles) < 1:
        print_error("There were no files found in the provided inputs that end with the accepted fastq extensions [blue].fq .fastq .fq.gz .fastq.gz[/blue]")
        sys.exit(1)

    # check if any names will be clashing
    bn_r = r"[\.\_](?:[RF])?(?:[12])?(?:\_00[1-9])*?$"
    uniqs = set()
    dupes = []
    for i in infiles:
        sans_ext = os.path.basename(re_ext.sub("", str(i)))
        if sans_ext in uniqs:
            dupes.append(sans_ext)
        else:
            uniqs.add(sans_ext)
    if dupes:
        print_error(Markdown("Identical filenames were detected, which will cause unexpected behavior and results.\n- files with identical names but different-cased extensions are treated as identical\n- files with the same name from different directories are also considered identical"))
        print_solution_with_culprits("Make sure all input files have unique names.", "Files with clashing names:")
        for i in dupes:
            click.echo(" ".join([j for j in infiles if i in j]), file = sys.stderr)
        sys.exit(1)

    n = len({re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in uniqs})
    # return the filenames and # of unique samplenames
    return infiles, n

def parse_alignment_inputs(inputs):
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
                    bam_infiles.append(Path(os.path.join(i, j)).resolve())
                elif re_bai.match(j):
                    bai_infiles.append(Path(os.path.join(i, j)).resolve())
        else:
            if re_bam.match(i):
                bam_infiles.append(Path(i).resolve())
            elif re_bai.match(i):
                bai_infiles.append(Path(i).resolve())
    if len(bam_infiles) < 1:
        print_error("There were no files found in the provided inputs that end with the [blue].bam[/blue] or [blue].sam[/blue] extensions.")
        sys.exit(1)
    re_ext = re.compile(r"\.(bam|sam)$", re.IGNORECASE)

    # check if any links will be clashing
    uniqs = set()
    dupes = []
    for i in bam_infiles:
        bn = os.path.basename(re_ext.sub("", str(i)))
        if bn in uniqs:
            dupes.append(bn)
        else:
            uniqs.add(bn)
    if dupes:
        print_error(Markdown("Identical filenames were detected, which will cause unexpected behavior and results.\n- files with identical names but different-cased extensions are treated as identical\n- files with the same name from different directories are also considered identical"))
        print_solution_with_culprits("Make sure all input files have unique names.", "Files with clashing names:")
        for i in dupes:
            click.echo(" ".join([j for j in bam_infiles if i in j]), file = sys.stderr)
        sys.exit(1)
    return bam_infiles, len(uniqs)

def biallelic_contigs(vcf, workdir):
    """Identify which contigs have at least 2 biallelic SNPs and write them to workdir/vcf.biallelic"""
    vbn = os.path.basename(vcf)
    os.makedirs(f"{workdir}/", exist_ok = True)
    valid = []
    vcfheader = subprocess.check_output(['bcftools', 'view', '-h', vcf]).decode().split('\n')
    header_contigs = [i.split(",")[0].replace("##contig=<ID=","") for i in vcfheader if i.startswith("##contig=")]
    if vcf.lower().endswith("bcf") and not os.path.exists(f"{vcf}.csi"):
        subprocess.run(f"bcftools index {vcf}".split())
    if vcf.lower().endswith("vcf.gz") and not os.path.exists(f"{vcf}.csi"):
        subprocess.run(f"bcftools index --tbi {vcf}".split())

    for contig in header_contigs:
        # Use bcftools to count the number of biallelic SNPs in the contig
        viewcmd = subprocess.Popen(['bcftools', 'view', '-r', contig, '-v', 'snps', '-m2', '-M2', '-c', '2', vcf], stdout=subprocess.PIPE)
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
    return f"{workdir}/{vbn}.biallelic", len(valid)