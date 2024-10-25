"""Harpy demultiplex workflows"""

import os
import sys
import yaml
from pathlib import Path
from rich import box
from rich.table import Table
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, snakemake_log
from ._validations import validate_demuxschema

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
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/demultiplex/")
@click.option('-s', '--schema', required = True, type=click.Path(exists=True, dir_okay=False, readable=True), help = 'Tab-delimited file of sample\<tab\>barcode')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Demultiplex", show_default=True,  help = 'Output directory name')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of a container')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False,  help = 'Setup the workflow and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('R1_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument('R2_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument('I1_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument('I2_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True))
def gen1(r1_fq, r2_fq, i1_fq, i2_fq, output_dir, schema, threads, snakemake, skip_reports, quiet, hpc, conda, setup_only):
    """
    Demultiplex Generation I haplotagged FASTQ files

    Use the R1, R2, I2, and I2 FASTQ files provided by the sequencing facility as inputs (in that exact order) provided after the options. 
    The `--schema` must be **tab** (or space) delimited, have **no header** (i.e. no column names), and be in the format of `sample`\<tab\>`barcode`,
    where `barcode` is the C- beadtag assigned to the sample (.e.g. `C01`, `C02`, etc.)
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/demultiplex_gen1.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake:
        command += snakemake

    validate_demuxschema(schema)
    os.makedirs(f"{workflowdir}", exist_ok=True)
    fetch_rule(workflowdir, "demultiplex_gen1.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "demultiplex_gen1")
    configs = {
        "workflow" : "demultiplex gen1",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "workflow_call" : command.rstrip(),
        "reports" : {
            "skip": skip_reports
        },
        "inputs" : {
            "demultiplex_schema" : Path(schema).resolve().as_posix(),
            "R1": Path(r1_fq).resolve().as_posix(),
            "R2": Path(r2_fq).resolve().as_posix(),
            "I1": Path(i1_fq).resolve().as_posix(),
            "I2": Path(i2_fq).resolve().as_posix()
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding= "utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))

    create_conda_recipes()
    if setup_only:
        sys.exit(0)
    
    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column("value", justify="left")
    start_text.add_row("Haplotag Type:", "Generation I")
    start_text.add_row("Demultiplex Schema:", schema)
    start_text.add_row("Output Folder:", output_dir + "/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    launch_snakemake(command, "demultiplex_gen1", start_text, output_dir, sm_log, quiet, "workflow/demux.gen1.summary")

demultiplex.add_command(gen1)
