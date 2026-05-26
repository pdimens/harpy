"""Harpy align workflows"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTQfile, FASTAfile
from harpy.common.cli_params import BwaParams, StrobeAlignParams, SnakemakeParams
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow
from harpy.validation.fasta import FASTA
from harpy.validation.fastq import FASTQ

@click.group()
@click.help_option('--help', hidden = True)
def align():
    """
    Align sequences to a reference genome

    Provide an additional subcommand `bwa` or `strobe` to get more information on using
    those aligners. Both have comparable performance, but `strobe` is typically faster.
    The aligners are not linked-read aware, but the workflows ensure linked-read information
    is carried over to the alignment records.
    """

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/align/bwa/")
@click.option('-w', '--depth-window', panel = "Parameters", default = 50000, show_default = True, type = click.IntRange(min = 50), help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', panel = "Parameters", type = BwaParams(), help = 'Additional bwa mem parameters, in quotes')
@click.option('-u', '--keep-unmapped', panel = "Parameters",  is_flag = True, default = False, help = 'Include unmapped sequences in output')
@click.option('-q', '--min-quality', panel = "Parameters", default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to output')
@click.option('-d', '--molecule-distance', panel = "Parameters", default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Distance cutoff for molecule assignment (bp)')
@click.option('-O', '--output', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Align/bwa", show_default=True,  help = 'Output directory name')
@click.option('-@', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-U','--unlinked', panel = "Parameters", is_flag = True, default = False, help = "Treat input data as not linked reads")
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('-C', '--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('-H', '--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('-Q', '--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('-T', '--no-temp', hidden = True, panel = "Workflow Options", is_flag = True, default = False, help = 'Don\'t delete temporary files')
@click.option('-N', '--setup', panel = "Workflow Options",  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('-R', '--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('-S', '--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.help_option('--help', hidden = True)
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def bwa(reference, inputs, output, depth_window, unlinked, threads, keep_unmapped, extra_params, min_quality, molecule_distance, snakemake, skip_reports, quiet, hpc, clean, container, no_temp, setup):
    """
    Run the BWA-MEM2 alignment workflow on a reference and one or more FASTQ inputs.
    
    Constructs and configures the `align_bwa` workflow, validates the reference and FASTQ inputs, attaches linked-read metadata, sets workflow parameters and notebooks, and initializes or sets up the workflow. Linked-read presence and type are auto-detected; set `unlinked` to disable detection. Setting `molecule_distance` greater than 0 enables alignment-distance–based barcode deconvolution for reporting only (barcodes are not modified). Input paths may be individual files or directories and can use shell wildcards.
    """
    workflow = Workflow("align_bwa", "align_bwa.smk", output, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake, no_temp)
    workflow.notebook_files = ["align_stats.ipynb", "align_lrstats.ipynb", "samtools_stats.ipynb"]
    workflow.conda = ["align", "qc"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc = not unlinked, quiet = quiet)
    fasta = FASTA(reference, quiet = quiet)

    workflow.linkedreads["type"] = fastq.lr_type
    workflow.linkedreads["standardized"] = {"BX" : fastq.bx_tag, "VX": fastq.vx_tag}
    workflow.notebooks["skip"] = skip_reports
    workflow.input(fasta.file, "reference")
    workflow.input(fastq.files, "fastq")
    workflow.param(fastq.illumina_old, "illumina-format-old")
    workflow.param(molecule_distance, "distance-threshold")
    workflow.param(min_quality, "min-map-quality")
    workflow.param(keep_unmapped, "keep-unmapped")
    workflow.param(depth_window, "depth-windowsize")
    if extra_params:
        workflow.param(extra_params, "extra")

    workflow.info = {
        "Samples": fastq.count,
        "Linked-Read Type": fastq.lr_type,
        "Reference": os.path.basename(reference),
        "Output Folder" : os.path.relpath(output) + "/"
    }

    workflow.initialize(setup)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/align/strobe/")
@click.option('-w', '--depth-window', panel = "Parameters", default = 50000, show_default = True, type = click.IntRange(min = 50), help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', panel = "Parameters", type = StrobeAlignParams(), help = 'Additional strobealign parameters, in quotes')
@click.option('-u', '--keep-unmapped', panel = "Parameters",  is_flag = True, default = False, help = 'Include unmapped sequences in output')
@click.option('-q', '--min-quality', panel = "Parameters", default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to output')
@click.option('-d', '--molecule-distance', panel = "Parameters", default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Distance cutoff for molecule assignment (bp)')
@click.option('-O', '--output', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Align/strobealign", show_default=True,  help = 'Output directory name')
@click.option('-@', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-U','--unlinked', panel = "Parameters", is_flag = True, default = False, help = "Treat input data as not linked reads")
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('-T', '--no-temp', hidden = True, panel = "Workflow Options", is_flag = True, default = False, help = 'Don\'t delete temporary files')
@click.option('-C', '--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('-H', '--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('-Q', '--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('-N', '--setup', panel = "Workflow Options",  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('-R', '--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('-S', '--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.help_option('--help', hidden = True)
@click.argument('reference', type=FASTAfile(), nargs = 1)
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def strobe(reference, inputs, output, unlinked, keep_unmapped, depth_window, threads, extra_params, min_quality, molecule_distance, snakemake, skip_reports, quiet, hpc, clean, container, no_temp, setup):
    """
    Run the Snakemake-based strobealign workflow to align input reads to a reference using strobealign.
    
    Detects linked-read metadata from inputs by default (can be disabled with `unlinked`) and accepts one or more FASTQ files or directories. If `molecule_distance` is greater than 0, alignment-distance based barcode deconvolution is enabled. Notebook generation, container, cleanup, and Snakemake execution are controlled via the corresponding arguments.
    
    Parameters:
        reference (str): Path to the reference FASTA file.
        inputs (list[str] | str): One or more input FASTQ files or directories (wildcards allowed).
        output (str): Output directory for workflow results.
        unlinked (bool): If True, ignore barcode detection and treat inputs as unlinked reads.
        keep_unmapped (bool): If True, retain unmapped reads in workflow outputs.
        depth_window (int): Window size used for depth calculations in reporting.
        threads (int): Number of threads to allocate for Snakemake.
        extra_params (dict | None): Additional workflow-specific parameters to pass through.
        min_quality (int): Minimum mapping quality threshold to record.
        molecule_distance (int): Distance threshold for molecule-based barcode deconvolution (0 disables).
        snakemake (str | None): Additional Snakemake options or configuration.
        skip_reports (bool): If True, skip generating notebook reports.
        quiet (bool): If True, reduce CLI verbosity.
        hpc (bool): If True, configure workflow for HPC execution.
        clean (bool): If True, remove intermediate files after run.
        container (str | None): Container image to use for workflow execution.
        no_temp (bool): If True, disable use of temporary directories by Snakemake.
        setup (bool): If True, perform workflow setup and exit without executing rules.
    """
    workflow = Workflow("align_strobe", "align_strobe.smk", output, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake, no_temp)
    workflow.notebook_files = ["align_stats.ipynb", "align_lrstats.ipynb", "samtools_stats.ipynb"]
    workflow.conda = ["align", "qc"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc= not unlinked, quiet = quiet)
    fasta = FASTA(reference, quiet = quiet)

    workflow.input(fasta.file, "reference")
    workflow.input(fastq.files,"fastq")
    workflow.linkedreads["type"] = fastq.lr_type
    workflow.linkedreads["standardized"] = {"BX" : fastq.bx_tag, "VX": fastq.vx_tag}
    workflow.param(fastq.illumina_old, "illumina-format-old")
    workflow.param(molecule_distance, "distance-threshold")
    workflow.param(min_quality, "min-map-quality")
    workflow.param(keep_unmapped, "keep-unmapped")
    workflow.param(depth_window, "depth-windowsize")
    if extra_params:
        workflow.param(extra_params, "extra")
    workflow.notebooks["skip"] = skip_reports

    workflow.info = {
        "Samples" : fastq.count,
        "Linked-Read Type" : fastq.lr_type,
        "Reference" : os.path.basename(reference),
        "Output Folder" : os.path.relpath(output) + "/"
    }

    workflow.initialize(setup)

align.add_command(bwa)
align.add_command(strobe)
