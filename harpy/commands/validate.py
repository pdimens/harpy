"""Harpy validation checks of FASTQ and BAM files for workflows"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTQfile, SAMfile
from harpy.common.cli_types_generic import PANEL_OPTIONS, SnakemakeParams
from harpy.common.printing import workflow_info
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow
from harpy.validation.sam import SAM
from harpy.validation.fastq import FASTQ

@click.group(options_metavar='', context_settings={"help_option_names" : []})
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
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ],
    "harpy validate fastq": [
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/validate/")
@click.rich_config(PANEL_OPTIONS)
@click.option_panel("Workflow Options", options = ["--help"],   panel_styles = {"border_style": "blue"})
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(1, 999, clamp = True), help = 'Number of threads to use')
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Validate/bam", show_default=True,  help = 'Output directory name')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only', panel = "Workflow Options",  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def bam(inputs, output_dir, threads, snakemake, quiet, hpc, container, setup_only):
    """
    Validate linked-read BAM file format

    Provide the input alignment (`.bam`) files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/betula*.bam`), or both.
    
    Validation checks if alignments have BX:Z: tags, that the barcodes are properly formatted for the
    auto-detected linked-read type, and that the filename matches the `@RG ID` tag. This **will not**
    fix your data, but it will report the number of records that feature errors to help you diagnose
    if file formatting will cause downstream issues. 
    """
    workflow = Workflow("validate_bam", "validate_bam.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["validate_bam.qmd"]
    workflow.conda = ["report"]

    ## checks and validations ##
    alignments = SAM(inputs, detect_bc=True, nonlinked_ok = False)

    workflow.config = {
        "workflow" : workflow.name,
        "linkedreads": {
            "type": alignments.lr_type
        },
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "inputs" : alignments.files
    }

    workflow.start_text = workflow_info(
        ("Alignment Files:", alignments.count),
        ("Barcode Type:", alignments.lr_type),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/validate/")
@click.rich_config(PANEL_OPTIONS)
@click.option_panel("Workflow Options", options = ["--help"],   panel_styles = {"border_style": "blue"})
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Validate/fastq", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(1, 999, clamp = True), help = 'Number of threads to use')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only', panel = "Workflow Options",  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def fastq(inputs, output_dir, threads, snakemake, quiet, hpc, container, setup_only):
    """
    Validate linked-read FASTQ file format

    Provide the input fastq files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/wombat*.fastq.gz`), or both.
    
    Validation checks if records in fastq files have properly formatted barcodes for the auto-detected
    linked-read type, and that any comments in the read headers conform to the SAM specification
    of `TAG:TYPE:VALUE`. This **will not** fix your data, but it will report the number of reads
    that feature errors to help you diagnose if file formatting will cause downstream issues. 
    """
    workflow = Workflow("validate_fastq", "validate_fastq.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["validate_fastq.qmd"]
    workflow.conda = ["report"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc=True, nonlinked_ok=False)

    ## setup workflow ##

    workflow.config = {
        "workflow" : workflow.name,
        "linkedreads": {
            "type": fastq.lr_type
        },
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "inputs" : fastq.files
    }

    workflow.start_text = workflow_info(
        ("FASTQ Files:", fastq.count),
        ("Linked-Read Type:", fastq.lr_type),
        ("Output Folder:", os.path.basename(output_dir) + "/"),

    )

    workflow.initialize(setup_only)

validate.add_command(bam)
validate.add_command(fastq)
