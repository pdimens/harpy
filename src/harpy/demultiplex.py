import rich_click as click
from .printfunctions import print_onstart
from .helperfunctions import fetch_rule
from .validations import validate_demuxschema
import sys
import os
import re


@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
def demultiplex():
    """
    Demultiplex haplotagged FASTQ files

    Check that you are using the correct haplotag method/technology, since the different
    barcoding approaches have very different demultiplexing strategies.

    **Haplotag Technologies**
    - `gen1`: the original haplotagging barcode strategy developed by Meier _et al._ (2021)
    """

docstring = {
    "harpy demultiplex gen1": [
        {
            "name": "Parameters",
            "options": ["--schema"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ]

}

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/demultiplex/")
@click.option('-s', '--schema', required = True, type=click.Path(exists=True, dir_okay=False), help = 'Tab-delimited file of sample\<tab\>barcode')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "Demultiplex", show_default=True, help = 'Name of output directory')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False,  help = 'Print the generated snakemake command and exit')
@click.argument('R1_FQ', required=True, type=click.Path(exists=True, dir_okay=False))
@click.argument('R2_FQ', required=True, type=click.Path(exists=True, dir_okay=False))
@click.argument('I1_FQ', required=True, type=click.Path(exists=True, dir_okay=False))
@click.argument('I2_FQ', required=True, type=click.Path(exists=True, dir_okay=False))
def gen1(r1_fq, r2_fq, i1_fq, i2_fq, output_dir, schema, threads, snakemake, skipreports, quiet, conda, print_only):
    """
    Demultiplex Generation I haplotagged FASTQ files

    Use the R1, R2, I2, and I2 FASTQ files provided by the sequencing facility as inputs (in that exact order) provided after the options. 
    The `--schema` must be **tab** (or space) delimited, have **no header** (i.e. no column names), and be in the format of `sample`\<tab\>`barcode`,
    where `barcode` is the C- beadtag assigned to the sample (.e.g. `C01`, `C02`, etc.)
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
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

    #check_demux_fastq(fastq_input)
    validate_demuxschema(schema)
    fetch_rule(workflowdir, "demultiplex.gen1.smk")

    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"R1: {r1_fq}\n")
        config.write(f"R2: {r2_fq}\n")
        config.write(f"I1: {i1_fq}\n")
        config.write(f"I2: {i2_fq}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samplefile: {schema}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    print_onstart(
        f"Haplotag Type: Generation I\nDemultiplex Schema: {schema}\nOutput Directory: {output_dir}",
        "demultiplex gen1"
    )
    return command