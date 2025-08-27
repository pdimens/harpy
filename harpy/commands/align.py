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

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def align():
    """
    Align sequences to a reference genome

    Provide an additional subcommand `bwa` or `strobe` to get more information on using
    those aligners. Both have comparable performance, but `strobe` is typically faster.
    The aligners are not linked-read aware, but the workflows
    """

module_docstring = {
    "harpy align": [
        {
            "name": "Commands",
            "commands": ["bwa", "strobe"],
            "panel_styles": {"border_style" : "blue"}
        }
    ]
}

docstring = {
    "harpy align bwa": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--keep-unmapped", "--molecule-distance", "--min-quality", "--unlinked"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--contigs", "--depth-window", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ],
    "harpy align strobe": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--keep-unmapped", "--molecule-distance", "--min-quality", "--read-length", "--unlinked"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--contigs", "--depth-window", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/align/bwa/")
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = click.IntRange(min = 50), help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', type = BwaParams(), help = 'Additional bwa mem parameters, in quotes')
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Include unmapped sequences in output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to output')
@click.option('-d', '--molecule-distance', default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Distance cutoff for molecule assignment (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Align/bwa", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-U','--unlinked', is_flag = True, default = False, help = "Treat input data as not linked reads")
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` unified progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def bwa(reference, inputs, output_dir, depth_window, unlinked, threads, keep_unmapped, extra_params, min_quality, molecule_distance, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Align sequences to reference genome using BWA MEM
 
    Provide the reference fasta followed by input fastq files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/echidna*.fastq.gz`), or both.
    
    BWA is a fast, robust, and reliable aligner that does not use barcodes when mapping.
    Presence and type of linked-read data is auto-detected, but can be deliberately ignored using `-U`.
    Setting `--molecule-distance` to `>0` activates alignment-distance based barcode deconvolution.
    """
    workflow = Workflow("align_bwa", "align_bwa.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["align_stats.qmd", "align_bxstats.qmd"]
    workflow.conda = ["align", "r", "qc"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc = not unlinked)
    fasta = FASTA(reference)
    if contigs:
        fasta.match_contigs(contigs) 

    workflow.config = {
        "workflow" : workflow.name,
        "alignment_quality" : min_quality,
        "keep_unmapped" : keep_unmapped,
        "depth_windowsize" : depth_window,
        "linkedreads": {
            "type" : fastq.lr_type == "none",
            "standardized": fastq.bx_tag,
            "distance_threshold" : molecule_distance,
        },
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "reference": fasta.file,
            "fastq": fastq.files
        }
    }

    workflow.start_text = workflow_info(
        ("Samples:", fastq.count),\
        ("Linked-Read Type:", fastq.lr_type),
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/align/strobe/")
@click.option('-w', '--depth-window', default = 50000, show_default = True, type = click.IntRange(min = 50), help = 'Interval size (in bp) for depth stats')
@click.option('-x', '--extra-params', type = StrobeAlignParams(), help = 'Additional strobealign parameters, in quotes')
@click.option('-u', '--keep-unmapped',  is_flag = True, default = False, help = 'Include unmapped sequences in output')
@click.option('-q', '--min-quality', default = 30, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality to output')
@click.option('-d', '--molecule-distance', default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Distance cutoff for molecule assignment (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Align/strobealign", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-U','--unlinked', is_flag = True, default = False, help = "Treat input data as not linked reads")
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` unified progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=FASTAfile(), nargs = 1)
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def strobe(reference, inputs, output_dir, unlinked, keep_unmapped, depth_window, threads, extra_params, min_quality, molecule_distance, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Align sequences to reference genome using strobealign
 
    Provide the reference fasta followed by the input fastq files and/or directories at the end of the command as individual
    files/folders, using shell wildcards (e.g. `data/echidna*.fastq.gz`), or both.
    
    strobealign is an ultra-fast aligner comparable to bwa for sequences >100bp and does 
    not use barcodes when mapping. Presence and type of linked-read data is auto-detected,
    but can be deliberately ignored using `-U`. Setting `--molecule-distance` to `>0` activates
    alignment-distance based barcode deconvolution.
    """
    workflow = Workflow("align_strobe", "align_strobe.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["align_stats.qmd", "align_bxstats.qmd"]
    workflow.conda = ["align", "r", "qc"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc= not unlinked)
    fasta = FASTA(reference)

    if contigs:
        fasta.match_contigs(contigs)

    workflow.config = {
        "workflow" : workflow.name,
        "alignment_quality" : min_quality,
        "keep_unmapped" : keep_unmapped,
        "depth_windowsize" : depth_window,
        "linkedreads": {
            "type" : fastq.lr_type == "none",
            "standardized": fastq.bx_tag,
            "distance_threshold" : molecule_distance,
        },
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "reference": fasta.file,
            "fastq": fastq.files
        }
    }

    workflow.start_text = workflow_info(
        ("Samples:", fastq.count),
        ("Linked-Read Type:", fastq.lr_type),
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

align.add_command(bwa)
align.add_command(strobe)
