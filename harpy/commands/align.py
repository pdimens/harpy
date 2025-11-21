"""Harpy align workflows"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTQfile, FASTAfile
from harpy.common.cli_types_generic import ContigList, SnakemakeParams
from harpy.common.cli_types_params import BwaParams, StrobeAlignParams
from harpy.common.printing import workflow_info
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow
from harpy.validation.fasta import FASTA
from harpy.validation.fastq import FASTQ

@click.group(context_settings={"help_option_names" : []})
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
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Align/bwa", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-U','--unlinked', panel = "Parameters", is_flag = True, default = False, help = "Treat input data as not linked reads")
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--contigs',  panel = "Workflow Options", type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup', panel = "Workflow Options",  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def bwa(reference, inputs, output_dir, depth_window, unlinked, threads, keep_unmapped, extra_params, min_quality, molecule_distance, snakemake, skip_reports, quiet, hpc, clean, container, contigs, setup):
    """
    Align sequences to reference genome using BWA MEM2
    
    Provide the reference fasta followed by input fastq files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/echidna*.fastq.gz`), or both.
    
    BWA is a fast, robust, and reliable aligner that does not use barcodes when mapping.
    Presence and type of linked-read data is auto-detected, but can be deliberately ignored using `-U`.
    Setting `--molecule-distance` to `>0` activates alignment-distance based barcode deconvolution.
    """
    workflow = Workflow("align_bwa", "align_bwa.smk", output_dir, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake)
    workflow.report_files = ["align_stats.qmd", "align_bxstats.qmd"]
    workflow.conda = ["align", "report", "qc"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc = not unlinked, quiet = quiet > 0)
    fasta = FASTA(reference, quiet = quiet > 0)
    if contigs:
        fasta.match_contigs(contigs) 

    workflow.linkedreads["type"] = fastq.lr_type
    workflow.linkedreads["standardized"] = fastq.bx_tag
    workflow.reports["skip"] = skip_reports
    workflow.reports["plot-contigs"] = contigs if contigs else "default"
    workflow.input(fasta.file, "reference")
    workflow.input(fastq.files, "fastq")
    workflow.param(molecule_distance, "distance-threshold")
    workflow.param(min_quality, "min-map-quality")
    workflow.param(keep_unmapped, "keep-unmapped")
    workflow.param(depth_window, "depth-windowsize")
    if extra_params:
        workflow.param(extra_params, "extra")

    workflow.start_text = workflow_info(
        ("Samples:", fastq.count),
        ("Linked-Read Type:", fastq.lr_type),
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.relpath(output_dir) + "/")
    )

    workflow.initialize(setup)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/align/strobe/")
@click.option('-w', '--depth-window', panel = "Parameters", default = 50000, show_default = True, type = click.IntRange(min = 50), help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', panel = "Parameters", type = StrobeAlignParams(), help = 'Additional strobealign parameters, in quotes')
@click.option('-u', '--keep-unmapped', panel = "Parameters",  is_flag = True, default = False, help = 'Include unmapped sequences in output')
@click.option('-q', '--min-quality', panel = "Parameters", default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to output')
@click.option('-d', '--molecule-distance', panel = "Parameters", default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Distance cutoff for molecule assignment (bp)')
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Align/strobealign", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-U','--unlinked', panel = "Parameters", is_flag = True, default = False, help = "Treat input data as not linked reads")
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--contigs', panel = "Workflow Options",  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup', panel = "Workflow Options",  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=FASTAfile(), nargs = 1)
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def strobe(reference, inputs, output_dir, unlinked, keep_unmapped, depth_window, threads, extra_params, min_quality, molecule_distance, snakemake, skip_reports, quiet, hpc, clean, container, contigs, setup):
    """
    Align sequences to reference genome using strobealign
 
    Provide the reference fasta followed by the input fastq files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/echidna*.fastq.gz`), or both.
    
    strobealign is an ultra-fast aligner comparable to bwa for sequences >100bp and does 
    not use barcodes when mapping. Presence and type of linked-read data is auto-detected,
    but can be deliberately ignored using `-U`. Setting `--molecule-distance` to `>0` activates
    alignment-distance based barcode deconvolution.
    """
    workflow = Workflow("align_strobe", "align_strobe.smk", output_dir, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake)
    workflow.report_files = ["align_stats.qmd", "align_bxstats.qmd"]
    workflow.conda = ["align", "report", "qc"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc= not unlinked, quiet= quiet > 0)
    fasta = FASTA(reference, quiet= quiet > 0)
    if contigs:
        fasta.match_contigs(contigs)

    workflow.input(fasta.file, "reference")
    workflow.input(fastq.files,"fastq")
    workflow.linkedreads["type"] = fastq.lr_type
    workflow.linkedreads["standardized"] = fastq.bx_tag
    workflow.param(molecule_distance, "distance-threshold")
    workflow.param(min_quality, "min-map-quality")
    workflow.param(keep_unmapped, "keep-unmapped")
    workflow.param(depth_window, "depth-windowsize")
    if extra_params:
        workflow.param(extra_params, "extra")
    workflow.reports["skip"] = skip_reports
    workflow.reports["plot-contigs"] = contigs if contigs else "default"

    workflow.start_text = workflow_info(
        ("Samples:", fastq.count),
        ("Linked-Read Type:", fastq.lr_type),
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.relpath(output_dir) + "/")
    )

    workflow.initialize(setup)

align.add_command(bwa)
align.add_command(strobe)
