import rich_click as click
from .helperfunctions import get_samples_from_fastq
import subprocess
import re
import os
import sys
import glob

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True)
@click.option('-d', '--directory', required = True, type=click.Path(exists=True), metavar = "Folder Path", help = 'Directory with raw sample sequences')
@click.option('-l', '--max-length', default = 150, show_default = True, type=int, metavar = "Integer", help = 'Maximum length to trim sequences down to')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional Fastp parameters, in quotes')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def qc(directory, max_length, extra_params, threads, snakemake, quiet):
    """
    Remove adapters and quality trim sequences
    """
    get_samples_from_fastq(directory)
    command = f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/qc.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    directory = directory.rstrip("/^")
    command.append(f"seq_directory={directory}")
    command.append(f"maxlen={max_length}")
    if extra_params is not None:
        command.append(f"extra={extra_params}")
    _module = subprocess.run(command)
    sys.exit(_module.returncode)