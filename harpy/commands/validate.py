"""Harpy validation checks of FASTQ and BAM files for workflows"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTQfile, SAMfile
from harpy.common.cli_types_generic import SnakemakeParams
from harpy.common.misc import container_ok
from harpy.validation.sam import SAM
from harpy.validation.fastq import FASTQ
from harpy.common.printing import workflow_info
from harpy.common.workflow import Workflow

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
            "options": ["--lr-type"],
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
            "options": ["--lr-type"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/validate/")
@click.option('-L', '--lr-type', type = click.Choice(['haplotagging', 'stlfr','tellseq'], case_sensitive=False), default = "haplotagging", show_default=True, help = "Linked read type\n[haplotagging, stlfr, tellseq]")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1, 999, clamp = True), help = 'Number of threads to use')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Validate/bam", show_default=True,  help = 'Output directory name')
@click.option('--quiet', show_default = True, default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` unified progress bar, `2` no output')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def bam(inputs, lr_type, output_dir, threads, snakemake, quiet, hpc, container, setup_only):
    """
    Run validity checks on haplotagged BAM files

    Provide the input alignment (`.bam`) files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/betula*.bam`), or both.
    
    It will check if alignments have BX:Z: tags, that the barcodes are properly formatted for the given
    `--lr-type` and that the filename matches the `@RG ID` tag.
    This **will not** fix your data, but it will report the number of records that feature errors  to help
    you diagnose if file formatting will cause downstream issues. 
    """
    workflow = Workflow("validate_bam", "validate_bam.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["validate_bam.qmd"]
    workflow.conda = ["r"]

    ## checks and validations ##
    alignments = SAM(inputs)

    workflow.config = {
        "workflow" : workflow.name,
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "linkedread_type": lr_type.lower(),
        "inputs" : alignments.files
    }

    workflow.start_text = workflow_info(
        ("Alignment Files:", alignments.count),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/validate/")
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Validate/fastq", show_default=True,  help = 'Output directory name')
@click.option('-L', '--lr-type', type = click.Choice(['haplotagging', 'stlfr','tellseq'], case_sensitive=False), default = "haplotagging", show_default=True, help = "Linked read type\n[haplotagging, stlfr, tellseq]")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1, 999, clamp = True), help = 'Number of threads to use')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` unified progress bar, `2` no output')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def fastq(inputs, output_dir, lr_type, threads, snakemake, quiet, hpc, container, setup_only):
    """
    Run validity checks on haplotagged FASTQ files.

    Provide the input fastq files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/wombat*.fastq.gz`), or both.
    
    It will check if fastq reads have properly formatted barcodes present in the location and format
    they are expected to be given `--lr-type` and that the comments in the read headers conform to the
    SAM specification of `TAG:TYPE:VALUE`. This **will not**
    fix your data, but it will report the number of reads that feature errors to help
    you diagnose if file formatting will cause downstream issues. 
    """
    workflow = Workflow("validate_fastq", "validate_fastq.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["validate_fastq.qmd"]
    workflow.conda = ["r"]

    ## checks and validations ##
    fastq = FASTQ(inputs)

    ## setup workflow ##

    workflow.config = {
        "workflow" : workflow.name,
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "linkedread_type": lr_type.lower(),
        "inputs" : fastq.files
    }

    workflow.start_text = workflow_info(
        ("FASTQ Files:", fastq.count),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

validate.add_command(bam)
validate.add_command(fastq)
