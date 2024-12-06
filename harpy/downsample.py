"""Downsample a fastq/bam file by barcodes"""

import os
import re
import sys
import pysam
import random
import subprocess
from pathlib import Path
import rich_click as click
from ._misc import harpy_pulsebar


docstring = {
    "harpy downsample": [
        {
            "name": "Parameters",
            "options": sorted(["--downsample", "--invalid", "--bx-tag", "--random-seed", "--prefix"]),
        },
        {
            "name": "Workflow Controls",
            "options": ["--quiet", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/downsample")
@click.option('-d', '--downsample', type = click.IntRange(min = 1), help = 'Downsampling amount')
@click.option('-i', '--invalid', default = "keep", show_default = True, type=click.Choice( ["keep","drop"]), help = "Strategy to handle invalid/missing barcodes")
@click.option('-b', '--bx-tag', type = str, default = "BX", show_default=True, help = "The header tag with the barcode")
@click.option('-p', '--prefix', type = click.Path(exists = False), default = "downsampled", show_default = True, help = 'Prefix for output file(s)')
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for sampling")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.argument('input', required=True, type=click.Path(exists=True, readable=True, dir_okay=False), nargs=-1)
def downsample(input, prefix, downsample, invalid, bx_tag, random_seed, threads, quiet):
    """
    Downsample data by barcode

    Downsamples FASTQ or BAM file(s) by barcode to keep all reads
    from `-d` barcodes.
    - one BAM file
    - two FASTQ files (R1 and R2 of a paired-end read set)
    
    Use `--invalid` to specify how to handle invalid barcodes:
    - `keep`: output all invalid/missing barcodes
    - `drop`: don't output any invalid/missing barcodes
    """
    # validate input files as either 1 bam or 2 fastq
    if len(bx_tag) != 2 or not bx_tag.isalnum():
        raise click.BadParameter(f'\'{bx_tag}\' is not a valid SAM tag. Tags for --bx-tag must be alphanumeric and exactly 2 characters, e.g. "BX"')
    if len(input) > 2:
        raise click.BadParameter('inputs must be 1 BAM file or 2 FASTQ files.')
    if len(input) == 1:
        if not input[0].lower().endswith(".bam"):
            raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.')
    if len(input) == 2:
        if input[0] == input[1]:
            raise click.BadParameter('the two input files cannot be identical')
        re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
        badfastq = []
        for i in input:
            if not re.search(re_ext, i):
                badfastq.append(i)
        if badfastq:
            raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.')

    output_bam = f"{prefix}.bam"
    bx_tag = bx_tag.upper()
    rng = random.Random(random_seed).sample if random_seed else random.Random().sample
    logfile = open(f"{prefix}.log", "w")
    sys.stderr.write = logfile.write

    with harpy_pulsebar(quiet, f"sorting by {bx_tag} tags") as progress:
        progress.add_task(f"sorting by {bx_tag} tags", total = None)
        if len(input) == 2:
            samtools_import = subprocess.Popen(f"samtools import -@ {threads} -T * {input[0]} {input[1]}".split(), stdout=subprocess.PIPE, stderr=logfile) 
            subprocess.run(f"samtools sort -@ {threads} -o bx_sorted.bam -t {bx_tag}".split(), stdin = samtools_import.stdout, stderr=logfile)
        else:
            subprocess.run(f"samtools sort -@ {threads} -o bx_sorted.bam -t {bx_tag} {input[0]}".split(), stderr=logfile)
    sorted_bam = Path("bx_sorted.bam").resolve().as_posix()
    invalid_pattern = re.compile(r'[AaBbCcDd]00')

    # read input file, get list of valid barcodes, subsample valid barcodes, read input again, output only sampled barcodes
    with (
        pysam.AlignmentFile(sorted_bam, "rb", check_sq=False) as infile,
        harpy_pulsebar(quiet, f"parsing {bx_tag} tags") as progress
    ):
        progress.add_task(f"parsing {bx_tag} tags", total = None)
        barcodes = set()
        for record in infile:
            try:
                barcode = record.get_tag(bx_tag)
                if isinstance(barcode, int):
                    pass # an int from an MI-tharype tag
                elif invalid_pattern.search(barcode):
                    continue
            except KeyError:
                continue
            barcodes.add(barcode)
        n_bc = len(barcodes)
        if downsample > n_bc:
            sys.stderr.write(f"The number of intended barcodes to downsample to ({downsample}) is greater than the number of barcodes in the input file ({n_bc}). Please choose a smaller number to downsample to.\n")
            sys.exit(1)
        barcodes = rng(sorted(barcodes), downsample)
    with (
        pysam.AlignmentFile(sorted_bam, "rb", check_sq=False) as infile,
        pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile,
        harpy_pulsebar(quiet, f"downsampling {bx_tag} tags") as progress
    ):
        progress.add_task(f"downsampling {bx_tag} tags", total = None)
        for record in infile:
            try:
                barcode = record.get_tag(bx_tag)
                if isinstance(barcode, int):
                    pass # an int from an MI-type tag
                elif invalid_pattern.search(barcode):
                    if invalid == "keep":
                        outfile.write(record)
                    continue
            except KeyError:
                if invalid == "keep":
                    outfile.write(record)
                continue
            if barcode in barcodes:
                outfile.write(record)

    # convert back
    os.remove(sorted_bam)
    if len(input) == 2:
        subprocess.run(f'samtools fastq -T * -c 6 -1 {prefix}.R1.fq.gz -2 {prefix}.R2.fq.gz {output_bam}'.split(), stderr=logfile)
        os.remove(output_bam)
    logfile.close()



# store records with the same barcode in an array
# then subsample according to the downsample fraction and write to output
# "within":
#        with (
#            pysam.AlignmentFile(sorted_bam, "rb", check_sq=False) as infile,
#            pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile,
#            harpy_pulsebar(quiet, f"downsampling within {bx_tag} tags") as progress
#        ):
#            progress.add_task(f"downsampling within {bx_tag} tags", total = None)
#            record_store_F = []
#            record_store_R = []
#            last_barcode = None
#            for record in infile:
#                try:
#                    barcode = record.get_tag(bx_tag)
#                    if isinstance(barcode, int):
#                        pass
#                    elif invalid_pattern.search(barcode):
#                        if invalid == "keep":
#                            outfile.write(record)
#                        elif invalid == "downsample":
#                            if rng.uniform(0, 1) <= downsample:
#                                outfile.write(record)    
#                        continue
#                except KeyError:
#                    if invalid == "keep":
#                        outfile.write(record)
#                    elif invalid == "downsample":
#                        if rng.uniform(0, 1) <= downsample:
#                            outfile.write(record)    
#                    continue
#                    
#                if last_barcode and barcode != last_barcode:
#                    # subsample records to file and reset record stores
#                    for i,j in zip(record_store_F, record_store_R):
#                        if rng.uniform(0, 1) <= downsample:
#                            if i:
#                                outfile.write(i)
#                            if j:
#                                outfile.write(j)
#                    record_store_F = []
#                    record_store_R = []
#                # otherwise proceed to append record to stores
#                if record.is_forward and record.is_paired:
#                    record_store_F.append(record)
#                elif record.is_reverse and record.is_paired:
#                    record_store_R.append(record)    
#                elif record.is_forward and not record.is_paired:
#                    record_store_F.append(record)
#                    record_store_R.append(None)
#                elif record.is_reverse and not record.is_paired:
#                    record_store_F.append(None)
#                    record_store_R.append(record)
#                # update what the last barcode is
#                last_barcode = barcode
