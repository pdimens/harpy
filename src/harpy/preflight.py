import rich_click as click
from .helperfunctions import fetch_file, generate_conda_deps, print_onstart
from .helperfunctions import parse_alignment_inputs, parse_fastq_inputs
import subprocess
import re
import os
import sys
import glob

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/preflight/")
#@click.option('-o', '--output', default = "current directory", show_default=True, type=click.Path(exists=False, writable=True, file_okay=False), metavar = "Folder Path", help = 'Output directory')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def fastq(input, threads, snakemake, quiet, print_only):
    """
    Run validity checks on haplotagged FASTQ files.

    For FASTQ sequence files, it will check that reads have `BX:Z:` tags, that haplotag
    barcodes are propery formatted (`AxxCxxBxxDxx`) and that the comments in the
    read headers conform to the SAM specification of `TAG:TYPE:VALUE`. This **will not**
    fix your data, but it will report the number of reads that feature errors to help
    you diagnose if file formatting will cause downstream issues. Provide the input fastq
    files and/or directories at the end of the command as either individual files/folders
    or using shell wildcards (e.g. `data/wombat*.fastq.gz`).
    """
    command = f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'Preflight/fastq/workflow/preflight-fastq.smk')
    command.append('--configfile')
    command.append(f"Preflight/fastq/workflow/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit()
    
    os.makedirs("Preflight/bam/workflow/", exist_ok= True)
    fetch_file("preflight-fastq.smk", f"Preflight/fastq/workflow/")
    fetch_file("PreflightFastq.Rmd", f"Preflight/fastq/workflow/report/")

    with open(f"Preflight/fastq/workflow/config.yml", "w") as config:
        config.write(f"seq_directory: Preflight/fastq/workflow/input\n")
        config.write(f"workflow_call: {call_SM}\n")

    sn = parse_fastq_inputs(input, "Preflight/fastq/workflow/input")

    print_onstart(
        f"Files: {len(sn)}\nOutput Directory: Preflight/fastq",
        "preflight fastq"
    )
    generate_conda_deps()
    _module = subprocess.run(command)
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/preflight/")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def bam(input, threads, snakemake, quiet, print_only):
    """
    Run validity checks on haplotagged BAM files

    Files must end in `.bam`. For BAM alignment files, it will check that alignments have BX:Z: tags, that haplotag
    barcodes are properly formatted (`AxxCxxBxxDxx`) and that the filename matches the `@RG ID` tag.
    This **will not** fix your data, but it will report the number of records that feature errors  to help
    you diagnose if file formatting will cause downstream issues. rovide the input fastq files and/or directories
    at the end of the command as either individual files/folders or using shell wildcards
    (e.g. `data/betula*.bam`).
    """
    command = f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'Preflight/bam/workflow/preflight-bam.smk')
    command.append('--configfile')
    command.append(f"Preflight/bam/workflow/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit()

    os.makedirs("Preflight/bam/workflow/", exist_ok= True)
    fetch_file("preflight-bam.smk", f"Preflight/bam/workflow/")
    fetch_file("PreflightBam.Rmd", f"Preflight/bam/workflow/report/")

    with open(f"Preflight/bam/workflow/config.yml", "w") as config:
        config.write(f"seq_directory: Preflight/bam/workflow/input\n")
        config.write(f"workflow_call: {call_SM}\n")

    sn = parse_alignment_inputs(input, "Preflight/bam/workflow/input")

    print_onstart(
        f"Samples: {len(sn)}\nOutput Directory: Preflight/bam",
        "preflight bam"
    )
    generate_conda_deps()
    _module = subprocess.run(command)
    sys.exit(_module.returncode)