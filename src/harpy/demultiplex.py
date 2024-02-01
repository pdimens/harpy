import rich_click as click
from .helperfunctions import fetch_file, generate_conda_deps, print_onstart
from .helperfunctions import validate_demuxschema, check_demux_fastq
import subprocess
import sys
import os
import re

@click.command(no_args_is_help = True)
@click.option('-f', '--file', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'The forward (or reverse) multiplexed FASTQ file')
@click.option('-b', '--samplesheet', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'Tab-delimited file of sample<tab>barcode')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
def gen1(file, samplesheet, threads, snakemake, skipreports, quiet, print_only):
    """
    Demultiplex Generation 1 haplotagged FASTQ files

    Use one of the four gzipped FASTQ files provided by the sequencer (I1, I2, R1, R2).fastq.gz for the `--file` argument, 
    Harpy will infer the other three. Note: the `--samplesheet` must be tab (or space) delimited and have no header (i.e. no column names).
    """
    check_demux_fastq(file)
    inprefix = re.sub(r"[\_\.][IR][12]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", "", os.path.basename(file))
    fetch_file("demultiplex.gen1.smk", f"Demultiplex/{inprefix}/workflow/")
    validate_demuxschema(samplesheet)

    command = f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'Demultiplex/{inprefix}/workflow/demultiplex.gen1.smk')
    command.append("--configfile")
    command.append(f"Demultiplex/{inprefix}/workflow/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]

    call_SM = " ".join(command)

    with open(f"Demultiplex/{inprefix}/workflow/config.yml", "w") as config:
        config.write(f"infile: {file}\n")
        config.write(f"samplefile: {samplesheet}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    if print_only:
        click.echo(call_SM)
    else:
        print_onstart(
            f"Initializing the [bold]harpy demultiplex gen1[/bold] workflow.\nInput Prefix: {inprefix}\nDemultiplex Schema: {samplesheet}"
        )
        generate_conda_deps()
        _module = subprocess.run(command)
        sys.exit(_module.returncode)