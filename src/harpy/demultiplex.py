"""Harpy demultiplex workflows"""

import os
import sys
import subprocess
from pathlib import Path
import rich_click as click
from .conda_deps import generate_conda_deps
from .printfunctions import print_onstart
from .helperfunctions import fetch_rule
from .validations import validate_demuxschema

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
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
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--skipreports", "--snakemake", "--threads", "--help"],
        },
    ]

}

@click.command(no_args_is_help = True, epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/demultiplex/")
@click.option('-s', '--schema', required = True, type=click.Path(exists=True, dir_okay=False, readable=True), help = 'Tab-delimited file of sample\<tab\>barcode')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Demultiplex", show_default=True,  help = 'Output directory name')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--config-only',  is_flag = True, hidden = True, default = False,  help = 'Create the config.yaml file and exit')
@click.argument('R1_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument('R2_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument('I1_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument('I2_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
def gen1(r1_fq, r2_fq, i1_fq, i2_fq, output_dir, schema, threads, snakemake, skipreports, quiet, hpc, conda, config_only):
    """
    Demultiplex Generation I haplotagged FASTQ files

    Use the R1, R2, I2, and I2 FASTQ files provided by the sequencing facility as inputs (in that exact order) provided after the options. 
    The `--schema` must be **tab** (or space) delimited, have **no header** (i.e. no column names), and be in the format of `sample`\<tab\>`barcode`,
    where `barcode` is the C- beadtag assigned to the sample (.e.g. `C01`, `C02`, etc.)
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/demultiplex_gen1.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake

    validate_demuxschema(schema)
    os.makedirs(f"{workflowdir}", exist_ok=True)
    fetch_rule(workflowdir, "demultiplex-gen1.smk")
        
    with open(f"{workflowdir}/config.yaml", "w", encoding= "utf-8") as config:
        config.write("workflow: demultiplex gen1\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"skip_reports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  demultiplex_schema: {Path(schema).resolve()}\n")
        config.write(f"  R1: {Path(r1_fq).resolve()}\n")
        config.write(f"  R2: {Path(r2_fq).resolve()}\n")
        config.write(f"  I1: {Path(i1_fq).resolve()}\n")
        config.write(f"  I2: {Path(i2_fq).resolve()}\n")
    if config_only:
        sys.exit(0)

    print_onstart(
        f"Haplotag Type: Generation I\nDemultiplex Schema: {schema}\nOutput Directory: {output_dir}",
        "demultiplex gen1"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)

demultiplex.add_command(gen1)
