import rich_click as click
from .printfunctions import print_onstart
from .helperfunctions import fetch_file, generate_conda_deps
from .validations import validate_demuxschema, check_demux_fastq
import subprocess
import sys
import os
import re

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/demultiplex/")
@click.option('-b', '--samplesheet', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'Tab-delimited file of sample\<tab\>barcode')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "Demultiplex", show_default=True, metavar = "String", help = 'Name of output directory')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('FASTQ_INPUT', required=True, type=click.Path(exists=True, dir_okay=False))
def gen1(fastq_input, output_dir, samplesheet, threads, snakemake, skipreports, quiet, print_only):
    """
    Demultiplex Generation 1 haplotagged FASTQ files

    Use one of the four gzipped FASTQ files provided by the sequencer (I1, I2, R1, R2).fastq.gz for the `FASTQ_INPUT` argument, 
    Harpy will infer the other three. Note: the `--samplesheet` must be tab (or space) delimited and have no header (i.e. no column names).
    """
    inprefix = re.sub(r"[\_\.][IR][12]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", "", os.path.basename(input))
    output_dir = output_dir.rstrip("/") + f"/{inprefix}"
    workflowdir = f"{output_dir}/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/demultiplex.gen1.smk')
    command.append("--configfile")
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

    check_demux_fastq(fastq_input)
    validate_demuxschema(samplesheet)
    fetch_file("demultiplex.gen1.smk", f"{workflowdir}/")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"infile: {fastq_input}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"infile_prefix: {inprefix}\n")
        config.write(f"samplefile: {samplesheet}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    print_onstart(
        f"Input Prefix: {inprefix}\nDemultiplex Schema: {samplesheet}\nOutput Directory: {output_dir}",
        "demultiplex gen1"
    )
    _module = subprocess.run(command)
    sys.exit(_module.returncode)