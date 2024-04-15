import rich_click as click
from .helperfunctions import fetch_rule, fetch_report, fetch_script, generate_conda_deps
from .printfunctions import print_onstart
from .fileparsers import parse_alignment_inputs, parse_fastq_inputs
import re
import os
import sys
import glob

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/preflight/")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "Preflight/fastq", show_default=True, metavar = "String", help = 'Name of output directory')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def fastq(input, output_dir, threads, snakemake, quiet, print_only):
    """
    Run validity checks on haplotagged FASTQ files.

    Provide the input fastq files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/wombat*.fastq.gz`), or both.
    
    It will check that fastq reads have `BX:Z:` tags, that haplotag
    barcodes are propery formatted (`AxxCxxBxxDxx`) and that the comments in the
    read headers conform to the SAM specification of `TAG:TYPE:VALUE`. This **will not**
    fix your data, but it will report the number of reads that feature errors to help
    you diagnose if file formatting will cause downstream issues. 
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/preflight-fastq.smk')
    command.append('--configfile')
    command.append(f"{workflowdir}/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit()
    
    os.makedirs(f"{workflowdir}/", exist_ok= True)
    sn = parse_fastq_inputs(input, f"{workflowdir}/input")
    fetch_rule(workflowdir, "preflight-fastq.smk")
    fetch_script(workflowdir, "checkFASTQ.py")
    fetch_report(workflowdir, "PreflightFastq.Rmd")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"workflow_call: {call_SM}\n")

    print_onstart(
        f"Files: {len(sn)}\nOutput Directory: {output_dir}/",
        "preflight fastq"
    )
    return command

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/preflight/")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "Preflight/bam", show_default=True, metavar = "String", help = 'Name of output directory')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def bam(input, output_dir, threads, snakemake, quiet, print_only):
    """
    Run validity checks on haplotagged BAM files

    Provide the input alignment (`.bam`) files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/betula*.bam`), or both.
    
    It will check that alignments have BX:Z: tags, that haplotag
    barcodes are properly formatted (`AxxCxxBxxDxx`) and that the filename matches the `@RG ID` tag.
    This **will not** fix your data, but it will report the number of records that feature errors  to help
    you diagnose if file formatting will cause downstream issues. 
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/preflight-bam.smk')
    command.append('--configfile')
    command.append(f"{workflowdir}/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit()

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    sn = parse_alignment_inputs(input, f"{workflowdir}/input")
    fetch_rule(workflowdir, "preflight-bam.smk")
    fetch_report(workflowdir, "PreflightBam.Rmd")
    fetch_script(workflowdir, "checkBAM.py")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"workflow_call: {call_SM}\n")

    print_onstart(
        f"Samples: {len(sn)}\nOutput Directory: {output_dir}/",
        "preflight bam"
    )
    return command