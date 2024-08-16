"""Harpy preflight-check workflows for FASTQ and BAM files"""

import os
import sys
from rich import box
from rich.table import Table
import rich_click as click
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_rule, fetch_report, snakemake_log, launch_snakemake
from .fileparsers import parse_alignment_inputs, parse_fastq_inputs

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def preflight():
    """
    File format checks for haplotag data

    This is useful to make sure your input files are formatted correctly for the processing pipeline 
    before you are surprised by errors hours into an analysis. Provide an additional command `fastq`
    or `bam` to see more information and options.
    """

docstring = {
    "harpy preflight bam": [
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
        },
    ],
    "harpy preflight fastq": [
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/preflight/")
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Preflight/fastq", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def fastq(inputs, output_dir, threads, snakemake, quiet, hpc, conda, setup_only):
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
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/preflight_fastq.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake is not None:
        command += snakemake

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    fqlist, n = parse_fastq_inputs(inputs)
    fetch_rule(workflowdir, "preflight_fastq.smk")
    fetch_report(workflowdir, "preflight_fastq.Rmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "preflight_fastq")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: preflight fastq\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        for i in fqlist:
            config.write(f"  - {i}\n")

    generate_conda_deps()
    if setup_only:
        sys.exit(0)

    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column("value", justify="left")
    start_text.add_row("FASTQ Files:", f"{n}")
    start_text.add_row("Output Folder:", output_dir + "/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", ""))
    launch_snakemake(command, "preflight_fastq", start_text, output_dir, sm_log, quiet)

@click.command(no_args_is_help = True, epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/preflight/")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Preflight/bam", show_default=True,  help = 'Output directory name')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def bam(inputs, output_dir, threads, snakemake, quiet, hpc, conda, setup_only):
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
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/preflight_bam.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake is not None:
        command += snakemake

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    bamlist, n = parse_alignment_inputs(inputs)
    fetch_rule(workflowdir, "preflight_bam.smk")
    fetch_report(workflowdir, "preflight_bam.Rmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "preflight_bam")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: preflight bam\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        for i in bamlist:
            config.write(f"  - {i}\n")

    generate_conda_deps()
    if setup_only:
        sys.exit(0)

    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column(header="value", justify="left")
    start_text.add_row("Alignment Files:", f"{n}")
    start_text.add_row("Output Folder:", output_dir + "/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", ""))
    launch_snakemake(command, "preflight_bam", start_text, output_dir, sm_log, quiet)

preflight.add_command(fastq)
preflight.add_command(bam)
