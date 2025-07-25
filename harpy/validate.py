"""Harpy validation checks of FASTQ and BAM files for workflows"""

import os
import sys
import yaml
import shutil
import rich_click as click
from .common.cli_types_generic import HPCProfile, SnakemakeParams
from .common.conda import create_conda_recipes
from .common.launch import launch_snakemake
from .common.misc import fetch_rule, fetch_report, instantiate_dir, setup_snakemake, write_workflow_config
from .common.parsers import parse_alignment_inputs, parse_fastq_inputs
from .common.printing import workflow_info

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def validate():
    """
    File format checks for linked-read data

    This is useful to make sure your input files are formatted correctly for the processing pipeline 
    before you are surprised by errors hours into an analysis. Provide an additional command `fastq`
    or `bam` to see more information and options.
    """

docstring = {
    "harpy validate bam": [
        {
            "name": "Parameters",
            "options": ["--platform"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ],
    "harpy validate fastq": [
        {
            "name": "Parameters",
            "options": ["--platform"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/validate/")
@click.option('-p', '--platform', type = click.Choice(['haplotagging', 'stlfr','tellseq'], case_sensitive=False), default = "haplotagging", show_default=True, help = "Linked read type\n[haplotagging, stlfr, tellseq]")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1, 999, clamp = True), help = 'Number of threads to use')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Validate/bam", show_default=True,  help = 'Output directory name')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def bam(inputs, platform, output_dir, threads, snakemake, quiet, hpc, container, setup_only):
    """
    Run validity checks on haplotagged BAM files

    Provide the input alignment (`.bam`) files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/betula*.bam`), or both.
    
    It will check if alignments have BX:Z: tags, that the barcodes are properly formatted for the given
    `--platform` and that the filename matches the `@RG ID` tag.
    This **will not** fix your data, but it will report the number of records that feature errors  to help
    you diagnose if file formatting will cause downstream issues. 
    """
    workflow = "validate_bam"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    bamlist, n = parse_alignment_inputs(inputs, "INPUTS")

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "validate_bam.smk")
    fetch_report(workflowdir, "validate_bam.qmd")

    conda_envs = ["r"]
    configs = {
        "workflow" : workflow,
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "platform": platform.lower(),
        "inputs" : bamlist
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Alignment Files:", n),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/validate.bam.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/validate/")
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Validate/fastq", show_default=True,  help = 'Output directory name')
@click.option('-p', '--platform', type = click.Choice(['haplotagging', 'stlfr','tellseq'], case_sensitive=False), default = "haplotagging", show_default=True, help = "Linked read type\n[haplotagging, stlfr, tellseq]")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1, 999, clamp = True), help = 'Number of threads to use')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable = True, resolve_path=True), nargs=-1)
def fastq(inputs, output_dir, platform, threads, snakemake, quiet, hpc, container, setup_only):
    """
    Run validity checks on haplotagged FASTQ files.

    Provide the input fastq files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/wombat*.fastq.gz`), or both.
    
    It will check if fastq reads have properly formatted barcodes present in the location and format
    they are expected to be given `--platform` and that the comments in the read headers conform to the
    SAM specification of `TAG:TYPE:VALUE`. This **will not**
    fix your data, but it will report the number of reads that feature errors to help
    you diagnose if file formatting will cause downstream issues. 
    """
    workflow = "validate_fastq"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    fqlist, n = parse_fastq_inputs(inputs, "INPUTS")

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "validate_fastq.smk")
    fetch_report(workflowdir, "validate_fastq.qmd")

    conda_envs = ["r"]
    configs = {
        "workflow" : workflow,
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "platform": platform.lower(),
        "inputs" : fqlist
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("FASTQ Files:", n),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/validate.fastq.summary")

validate.add_command(bam)
validate.add_command(fastq)
