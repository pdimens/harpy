"""Module for workflow file parsers"""

import sys
import os
import re
import glob
from pathlib import Path
import rich_click as click
from .printfunctions import print_error, print_solution_with_culprits, print_solution

def getnames(directory, ext):
    """Find all files in 'directory' that end with 'ext'"""
    samplenames = set([i.split(ext)[0] for i in os.listdir(directory) if i.endswith(ext)])
    if len(samplenames) < 1:
        print_error(f"No sample files ending with [bold]{ext}[/bold] found in [bold]{directory}[/bold].")
        sys.exit(1)
    return samplenames

def parse_fastq_inputs(inputs, outdir):
    """
    Parse the command line input FASTQ arguments to generate a clean list of input files
    and create symlinks of those files to a target destination folder.
    """
    infiles = []
    outfiles = []
    for i in inputs:
        if os.path.isdir(i):
            for j in os.listdir(i):
                if j.lower().endswith("gz") or j.lower().endswith("fastq") or j.lower().endswith("fq"):
                    infiles.append(os.path.join(i, j))
        else:
            if i.lower().endswith("gz") or i.lower().endswith("fastq") or i.lower().endswith("fq"):
                infiles.append(i)
    re_fq = re.compile(r"\.(fq|fastq)", re.IGNORECASE)
    re_gz = re.compile(r"\.gz$", re.IGNORECASE)
    re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
    if len(infiles) < 1:
        print_error("There were no files found in the provided inputs that end with the accepted fastq extensions [blue].fq .fastq .fq.gz .fastq.gz[/blue]")
        sys.exit(1)
    for i in infiles:
        destination = os.path.join(outdir,os.path.basename(i))
        # clean up extensions for consistency
        clean_destination = re_fq.sub(".fq", destination)
        clean_destination = re_gz.sub(".gz", clean_destination)
        outfiles.append(clean_destination)
    # check if any links will be clashing
    uniqs = set()
    dupes = [os.path.basename(re_ext.sub("", i)) for i in outfiles if i in uniqs or uniqs.add(i)]
    if dupes:
        print_error("Identical filenames were detected, which will cause unexpected behavior and results. Note that files with identical names but different-cased extensions are treated as identical.")
        print_solution_with_culprits("Make sure all input files have unique names.", "Files with clashing names:")
        for i in dupes:
            click.echo(" ".join([j for j in infiles if i in j]), file = sys.stderr)
        sys.exit(1)

    Path(outdir).mkdir(parents=True, exist_ok=True)
    for (i,o) in zip(infiles, outfiles):
        Path(o).unlink(missing_ok=True)
        Path(o).symlink_to(Path(i).resolve())
    return infiles

def parse_alignment_inputs(inputs, outdir):
    """
    Parse the command line input sam/bam arguments to generate a clean list of input files
    and create symlinks of those files to a target destination folder.
    """
    bam_infiles = []
    bai_infiles = []
    bam_outfiles = []
    bai_outfiles = []
    for i in inputs:
        if os.path.isdir(i):
            for j in os.listdir(i):
                if j.lower().endswith("bam") or j.lower().endswith("sam"):
                    bam_infiles.append(os.path.join(i, j))
                elif j.lower().endswith("bai"):
                    bai_infiles.append(os.path.join(i, j))
        else:
            if i.lower().endswith("bam") or i.lower().endswith("sam"):
                bam_infiles.append(i)
            elif i.lower().endswith("bai"):
                bai_infiles.append(i)
    if len(bam_infiles) < 1:
        print_error("There were no files found in the provided inputs that end with the [blue].bam[/blue] extension.")
        sys.exit(1)
    re_bam = re.compile(r"\.bam$", re.IGNORECASE)
    re_sam = re.compile(r"\.sam$", re.IGNORECASE)
    re_ext = re.compile(r"\.(bam|sam)$", re.IGNORECASE)
    for i in bam_infiles:
        destination = os.path.join(outdir,os.path.basename(i))
        # clean up extensions for consistency
        clean_destination = re_bam.sub(".bam", destination)
        clean_destination = re_sam.sub(".sam", clean_destination)
        bam_outfiles.append(clean_destination)
    # check if any links will be clashing
    uniqs = set()
    dupes = []
    for i in bam_infiles:
        bn = os.path.basename(re_ext.sub("", i))
        if bn in uniqs:
            dupes.append(bn)
        else:
            uniqs.add(bn)
    if dupes:
        print_error("Identical filenames were detected, which will cause unexpected behavior and results. Note that files with identical names but different-cased extensions are treated as identical.")
        print_solution_with_culprits("Make sure all input files have unique names.", "Files with clashing names:")
        for i in dupes:
            click.echo(" ".join([j for j in bam_infiles if i in j]), file = sys.stderr)
        sys.exit(1)
    for i in bai_infiles:
        destination = os.path.join(outdir,os.path.basename(i))
        # clean up extensions for consistency
        clean_destination = re_bam.sub(".bam.bai", destination)
        bai_outfiles.append(clean_destination)
    Path(outdir).mkdir(parents=True, exist_ok=True)
    for (i,o) in zip(bam_infiles, bam_outfiles):
        Path(o).unlink(missing_ok=True)
        Path(o).symlink_to(Path(i).resolve())
    for (i,o) in zip(bai_infiles, bai_outfiles):
        Path(o).unlink(missing_ok=True)
        Path(o).symlink_to(Path(i).resolve())
    return bam_infiles

def get_samples_from_fastq(directory):
    """Identify the sample names from a directory containing FASTQ files"""
    full_flist = [i for i in glob.iglob(f"{directory}/*") if not os.path.isdir(i)]
    r = re.compile(r".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    full_fqlist = list(filter(r.match, full_flist))
    fqlist = [os.path.basename(i) for i in full_fqlist]
    bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
    if len(fqlist) == 0:
        print_error(f"No fastq files with acceptable names found in [bold]{directory}[/bold]")
        print_solution("Check that the file endings conform to [green].[/green][[green]F[/green][dim]|[/dim][green]R1[/green]][green].[/green][[green]fastq[/green][dim]|[/dim][green]fq[/green]][green].gz[/green]\nRead the documentation for details: https://pdimens.github.io/harpy/haplotagdata/#naming-conventions")
        sys.exit(1)

    return set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])