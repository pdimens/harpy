"""Separate barcodes by unique molecule"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTQfile
from harpy.common.cli_params import SnakemakeParams
from harpy.validation.fastq import FASTQ
from harpy.common.system_ops import container_ok, is_arm
from harpy.common.workflow import Workflow

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/deconvolve")
@click.option('-k', '--kmer-length', panel = "Parameters", default = 21, show_default = True, type=click.IntRange(min = 1), help = 'Size of kmers')
@click.option('-w', '--window-size', panel = "Parameters", default = 40, show_default = True, type=click.IntRange(min = 3), help = 'Size of window guaranteed to contain at least one kmer')
@click.option('-d', '--density', panel = "Parameters", default = 3, show_default = True, type = click.IntRange(min = 1), help = 'On average, 1/2^d kmers are indexed')
@click.option('-a', '--dropout', panel = "Parameters", default = 0, show_default = True, type = click.IntRange(min = 0), help = 'Minimum cloud size to deconvolve')
@click.option('-@', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(1,999, clamp = True), help = 'Number of threads to use')
@click.option('-O', '--output', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Deconvolve", show_default=True,  help = 'Output directory name')
@click.option('-T', '--no-temp', hidden = True, panel = "Workflow Options", is_flag = True, default = False, help = 'Don\'t delete temporary files')
@click.option('-C', '--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('-H', '--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('-Q', '--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('-N', '--setup', panel = "Workflow Options",  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('-S', '--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.help_option('--help', hidden = True)
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def deconvolve(inputs, output, kmer_length, window_size, density, dropout, threads, snakemake, quiet, hpc, clean, container, setup, no_temp):
    """
    Resolve barcode sharing in unrelated molecules

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    
    The term "cloud" refers to the collection of all sequences that feature the same barcode. By default,
    `dropout` is set to `0`, meaning it will consider all barcodes, even clouds with singleton.
    """
    is_arm(allowed = False)
    workflow = Workflow("deconvolve", "deconvolve.smk", output, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake, no_temp)
    workflow.conda = ["qc"]

    ## checks and validations ##
    fastq = FASTQ(inputs, quiet= quiet > 0)
    
    workflow.input(fastq.files)
    workflow.param(kmer_length, "kmer-length")       
    workflow.param(window_size, "window-size")
    workflow.param(density, "density")
    workflow.param(dropout, "dropout")

    workflow.info = {
        "Samples:" : fastq.count,
        "Output Folder:" : os.path.relpath(output) + "/"
    }
    
    workflow.initialize(setup)
