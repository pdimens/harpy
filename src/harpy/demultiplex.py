import rich_click as click
from .helperfunctions import validate_demuxschema, check_demux_fastq
import subprocess
import sys
import os
import re

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True)
@click.option('-f', '--file', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'The forward (or reverse) multiplexed FASTQ file')
@click.option('-b', '--samplesheet', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'Tab-delimited file of sample<tab>barcode')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
def gen1(file, samplesheet, threads, snakemake, quiet):
    """
    Demultiplex Generation 1 haplotagged FASTQ files

    Use one of the four gzipped FASTQ files provided by the sequencer (I1, I2, R1, R2).fastq.gz for the `--file` argument, 
    Harpy will infer the other three. Note: the `--samplesheet` must be tab (or space) delimited and have no header (i.e. no column names).
    """
    check_demux_fastq(file)
    validate_demuxschema(samplesheet)

    command = f'snakemake --rerun-incomplete --nolock --cores {threads} --directory . --snakefile {harpypath}/demultiplex.gen1.smk'.split()
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    if quiet:
        command.append("--quiet")
        command.append("all")
    command.append('--config')
    command.append(f"infile={file}")
    command.append(f"samplefile={samplesheet}")
    _module = subprocess.run(command)
    sys.exit(_module.returncode)