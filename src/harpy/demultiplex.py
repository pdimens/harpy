"""Harpy demultiplex workflows"""

import os
import sys
import shutil
import subprocess
from pathlib import Path
import rich_click as click
from rich.markdown import Markdown
from .conda_deps import generate_conda_deps
from .printfunctions import print_onstart, print_notice
from .helperfunctions import fetch_rule, symlink
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
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ]

}

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/demultiplex/")
@click.option('-s', '--schema', required = True, type=click.Path(exists=True, dir_okay=False), help = 'Tab-delimited file of sample\<tab\>barcode')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "Demultiplex", show_default=True, help = 'Name of output directory')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False,  help = 'Print the generated snakemake command and exit')
@click.argument('R1_FQ', required=True, type=click.Path(exists=True, dir_okay=False))
@click.argument('R2_FQ', required=True, type=click.Path(exists=True, dir_okay=False))
@click.argument('I1_FQ', required=True, type=click.Path(exists=True, dir_okay=False))
@click.argument('I2_FQ', required=True, type=click.Path(exists=True, dir_okay=False))
def gen1(r1_fq, r2_fq, i1_fq, i2_fq, output_dir, schema, threads, snakemake, skipreports, quiet, hpc, conda, print_only):
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
    command += f"--snakefile {workflowdir}/demultiplex.gen1.smk "
    command += f"--configfile {workflowdir}/config.yml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake
    if print_only:
        click.echo(command)
        sys.exit(0)

    validate_demuxschema(schema)
    fetch_rule(workflowdir, "demultiplex.gen1.smk")
    os.makedirs(f"{workflowdir}/input", exist_ok=True)

    # copy or symlink
    Path(f"{workflowdir}/input/DATA_R1_001.fastq.gz").unlink(missing_ok=True)
    Path(f"{workflowdir}/input/DATA_R2_001.fastq.gz").unlink(missing_ok=True)
    Path(f"{workflowdir}/input/DATA_I1_001.fastq.gz").unlink(missing_ok=True)
    Path(f"{workflowdir}/input/DATA_I2_001.fastq.gz").unlink(missing_ok=True)
    if hpc:
        print_notice(Markdown(f"Copying input files into `{workflowdir}/input/`, which will delay starting the workflow."))
        shutil.copy(r1_fq, f"{workflowdir}/input/DATA_R1_001.fastq.gz")
        shutil.copy(r2_fq, f"{workflowdir}/input/DATA_R2_001.fastq.gz")
        shutil.copy(i1_fq, f"{workflowdir}/input/DATA_I1_001.fastq.gz")
        shutil.copy(i2_fq, f"{workflowdir}/input/DATA_I2_001.fastq.gz")
    else:
        symlink(r1_fq, f"{workflowdir}/input/DATA_R1_001.fastq.gz")
        symlink(r2_fq, f"{workflowdir}/input/DATA_R2_001.fastq.gz")
        symlink(i1_fq, f"{workflowdir}/input/DATA_I1_001.fastq.gz")
        symlink(i2_fq, f"{workflowdir}/input/DATA_I2_001.fastq.gz")
        
    with open(f"{workflowdir}/config.yml", "w", encoding= "utf-8") as config:
        config.write(f"R1: {workflowdir}/input/DATA_R1_001.fastq.gz\n")
        config.write(f"R2: {workflowdir}/input/DATA_R2_001.fastq.gz\n")
        config.write(f"I1: {workflowdir}/input/DATA_I1_001.fastq.gz\n")
        config.write(f"I2: {workflowdir}/input/DATA_I2_001.fastq.gz\n")
        #config.write(f"hpc: {bool(hpc)}")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samplefile: {schema}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")

    print_onstart(
        f"Haplotag Type: Generation I\nDemultiplex Schema: {schema}\nOutput Directory: {output_dir}",
        "demultiplex gen1"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)

demultiplex.add_command(gen1)
