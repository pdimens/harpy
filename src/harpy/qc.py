import rich_click as click
from .helperfunctions import generate_conda_deps, get_samples_from_fastq, fetch_file, print_onstart
import subprocess
import re
import os
import sys
import glob

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/qc")
@click.option('-d', '--directory', required = True, type=click.Path(exists=True, file_okay=False), metavar = "Folder Path", help = 'Directory with raw sample sequences')
@click.option('-l', '--max-length', default = 150, show_default = True, type=int, metavar = "Integer", help = 'Maximum length to trim sequences down to')
@click.option('-a', '--ignore-adapters', is_flag = True, show_default = False, default = False, metavar = "Toggle", help = 'Skip adapter trimming')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional Fastp parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def qc(directory, max_length, ignore_adapters, extra_params, threads, snakemake, skipreports, quiet, print_only):
    """
    Remove adapters and quality trim sequences

    By default, adapters will be automatically detected and removed (can be disabled).
    The input reads will be quality trimmed using:
    - a sliding window from front to tail
    - minimum 15bp length filter
    - poly-G tail removal
    """
    fetch_file("qc.smk", "QC/workflow/")
    fetch_file("BxCount.Rmd", "QC/workflow/report/")
    directory = directory.rstrip("/^")
    sn = get_samples_from_fastq(directory)

    command = f'snakemake --rerun-incomplete --nolock  --use-conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append('QC/workflow/qc.smk')
    command.append('--configfile')
    command.append('QC/workflow/config.yml')
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]

    call_SM = " ".join(command)

    with open("QC/workflow/config.yml", "w") as config:
        config.write(f"seq_directory: {directory}\n")
        config.write(f"adapters: {ignore_adapters}\n")
        config.write(f"maxlen: {max_length}\n")
        if extra_params is not None:
            command.append(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    if print_only:
        click.echo(call_SM)
    else:
        print_onstart(
            f"Input Directory: {directory}\nSamples: {len(sn)}",
            "qc"
        )
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)