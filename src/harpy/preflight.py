import rich_click as click
from .helperfunctions import fetch_file, generate_conda_deps, getnames, get_samples_from_fastq, print_onstart
import subprocess
import re
import os
import sys
import glob

@click.command(no_args_is_help = True)
@click.option('-d', '--directory', required = True, type=click.Path(exists=True, file_okay=False), metavar = "Folder Path", help = 'Directory with FASTQ files')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def fastq(directory, threads, snakemake, quiet, print_only):
    """
    Run validity checks on haplotagged FASTQ files.

    For FASTQ sequence files, it will check that reads have `BX:Z:` tags, that haplotag
    barcodes are propery formatted (`AxxCxxBxxDxx`) and that the comments in the
    read headers conform to the SAM specification of `TAG:TYPE:VALUE`. This **will not**
    fix your data, but it will report the number of reads that feature errors to help
    you diagnose if file formatting will cause downstream issues.
    """
    fetch_file("preflight-fastq.smk", f"{directory}/Preflight/workflow/")
    fetch_file("PreflightFastq.Rmd", f"{directory}/Preflight/workflow/report/")
    sn = get_samples_from_fastq(directory)
    directory = directory.rstrip("/^")
    command = f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{directory}/Preflight/workflow/preflight-fastq.smk')
    command.append('--configfile')
    command.append(f"{directory}/Preflight/workflow/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]

    call_SM = " ".join(command)

    with open(f"{directory}/Preflight/workflow/config.yml", "w") as config:
        config.write(f"seq_directory: {directory}\n")
        config.write(f"workflow_call: {call_SM}\n")

    if print_only:
        click.echo(call_SM)
    else:
        print_onstart(
            f"Initializing the [bold]harpy preflight fastq[/bold] workflow.\nInput Directory: {directory}\nSamples: {len(sn)}"
        )
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)

@click.command(no_args_is_help = True)
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with FASTQ files')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def bam(directory, threads, snakemake, quiet, print_only):
    """
    Run validity checks on haplotagged BAM files

    Files must end in `.bam` (lowercase). For BAM alignment files, it will check that alignments have BX:Z: tags, that haplotag
    barcodes are properly formatted (`AxxCxxBxxDxx`) and that the filename matches the `@RG ID` tag.
    This **will not** fix your data, but it will report the number of records that feature errors  to help
    you diagnose if file formatting will cause downstream issues.
    """
    fetch_file("preflight-bam.smk", f"{directory}/Preflight/workflow/")
    fetch_file("PreflightBam.Rmd", f"{directory}/Preflight/workflow/report/")
    directory = directory.rstrip("/^")
    flist = getnames(directory, ".bam")
    command = f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{directory}/Preflight/workflow/preflight-bam.smk')
    command.append('--configfile')
    command.append(f"{directory}/Preflight/workflow/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    
    call_SM = " ".join(command)

    with open(f"{directory}/Preflight/workflow/config.yml", "w") as config:
        config.write(f"seq_directory: {directory}\n")
        config.write(f"workflow_call: {call_SM}\n")

    if print_only:
        click.echo(call_SM)
    else:
        print_onstart(
            f"Initializing the [bold]harpy preflight bam[/bold] workflow.\nInput Directory: {directory}\nFiles: {len(flist)}"
        )
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)