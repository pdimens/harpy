import rich_click as click
from .helperfunctions import generate_conda_deps, get_samples_from_fastq, fetch_file, print_onstart, parse_fastq_inputs
import subprocess
import re
import os
import sys
import glob

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/qc")
@click.option('-l', '--max-length', default = 150, show_default = True, type=int, metavar = "Integer", help = 'Maximum length to trim sequences down to')
@click.option('-a', '--ignore-adapters', is_flag = True, show_default = False, default = False, metavar = "Toggle", help = 'Skip adapter trimming')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional Fastp parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def qc(input, max_length, ignore_adapters, extra_params, threads, snakemake, skipreports, quiet, print_only):
    """
    Remove adapters and quality trim sequences

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    
    By default, adapters will be automatically detected and removed (can be disabled).
    The input reads will be quality trimmed using:
    - a sliding window from front to tail
    - minimum 15bp length filter
    - poly-G tail removal
    """
    workflowdir = "QC/workflow"
    command = f'snakemake --rerun-incomplete --nolock  --use-conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/qc.smk')
    command.append('--configfile')
    command.append(f'{workflowdir}/config.yml')
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit(0)

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    _ = parse_fastq_inputs(input, f"{workflowdir}/input")
    sn = get_samples_from_fastq(f"{workflowdir}/input")
    fetch_file("qc.smk", f"{workflowdir}/")
    fetch_file("BxCount.Rmd", f"{workflowdir}/report/")
    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"adapters: {ignore_adapters}\n")
        config.write(f"maxlen: {max_length}\n")
        if extra_params is not None:
            command.append(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    print_onstart(
        f"Samples: {len(sn)}\nOutput Directory: QC/",
        "qc"
    )
    _module = subprocess.run(command)
    sys.exit(success)