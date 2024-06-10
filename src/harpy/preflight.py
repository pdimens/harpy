"""Harpy preflight-check workflows for FASTQ and BAM files"""

import os
import sys
import subprocess
import rich_click as click
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_rule, fetch_report, fetch_script
from .printfunctions import print_onstart
from .fileparsers import parse_alignment_inputs, parse_fastq_inputs

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def preflight():
    """
    Run file format checks on haplotag data

    This is useful to make sure your input files are formatted correctly for the processing pipeline 
    before you are surprised by errors hours into an analysis. Provide an additional command `fastq`
    or `bam` to see more information and options.
    """

docstring = {
    "harpy preflight bam": [
        {
            "name": "Options",
            "options": ["--output-dir", "--threads", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy preflight fastq": [
        {
            "name": "Options",
            "options": ["--output-dir", "--threads", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/preflight/")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "Preflight/fastq", show_default=True, help = 'Name of output directory')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def fastq(inputs, output_dir, threads, snakemake, quiet, hpc, conda, print_only):
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
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/preflight-fastq.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    if print_only:
        click.echo(command)
        sys.exit()
    
    os.makedirs(f"{workflowdir}/", exist_ok= True)
    sn = parse_fastq_inputs(inputs, f"{workflowdir}/input", hpc)
    fetch_rule(workflowdir, "preflight-fastq.smk")
    fetch_script(workflowdir, "checkFASTQ.py")
    fetch_report(workflowdir, "PreflightFastq.Rmd")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: preflight fastq\n")
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"workflow_call: {command}\n")

    print_onstart(
        f"Files: {len(sn)}\nOutput Directory: {output_dir}/",
        "preflight fastq"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/preflight/")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "Preflight/bam", show_default=True, help = 'Name of output directory')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def bam(inputs, output_dir, threads, snakemake, quiet, hpc, conda, print_only):
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
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/preflight-bam.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    if print_only:
        click.echo(command)
        sys.exit()

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    sn = parse_alignment_inputs(inputs, f"{workflowdir}/input", hpc)
    fetch_rule(workflowdir, "preflight-bam.smk")
    fetch_report(workflowdir, "PreflightBam.Rmd")
    fetch_script(workflowdir, "checkBAM.py")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: preflight bam\n")
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"workflow_call: {command}\n")

    print_onstart(
        f"Samples: {len(sn)}\nOutput Directory: {output_dir}/",
        "preflight bam"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)

preflight.add_command(fastq)
preflight.add_command(bam)
