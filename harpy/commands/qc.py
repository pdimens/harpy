"""Harpy sequence adapter trimming and quality control"""

import os

import rich_click as click

from harpy.common.cli_filetypes import FASTQfile, HPCProfile, QCAdapters
from harpy.common.cli_params import FastpParams, MultiInt, SnakemakeParams
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow
from harpy.validation.fastq import FASTQ

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/qc")
@click.option('-d', '--deduplicate', panel = "Parameters", is_flag = True, default = False, help = 'Identify and remove PCR duplicates')
@click.option('-x', '--extra-params', panel = "Parameters", type = FastpParams(), help = 'Additional Fastp parameters, in quotes')
@click.option('-l', '--length', panel = "Parameters", default = "30,150", type=MultiInt(2, minimum=30), show_default = True, help = 'Minimum,Maximum read lengths')
@click.option('-a', '--trim-adapters', panel = "Parameters", type = QCAdapters(), help = 'Find and remove adapter sequences')
@click.option('-U', '--unlinked', panel = "Parameters", is_flag = True, default = False, help = "Treat input data as not linked reads")
@click.option('-O', '--output', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "QC", show_default=True,  help = 'Output directory name')
@click.option('-@', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-T', '--no-temp', hidden = True, panel = "Workflow Options", is_flag = True, default = False, help = 'Don\'t delete temporary files')
@click.option('-C', '--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('-N', '--setup', panel = "Workflow Options",  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('-H', '--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('-Q', '--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('-R', '--skip-reports', panel = "Workflow Options",  is_flag = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('-S', '--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.help_option('--help', hidden = True)
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def qc(inputs, output, unlinked,length, trim_adapters, deduplicate, extra_params, threads, snakemake, skip_reports, quiet, hpc, clean, container, setup, no_temp):
    """
    FASTQ adapter removal, quality filtering, etc.

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    Linked-read presence and type is auto-detected, but you may use `-U` to disable the parts
    of the workflow specific to linked-read data. Adapter removal is optional and accepts `auto`
    for automatic detection, a nucleotide sequence, or a FASTA file of adapters to remove.

    **Mandatory QC**
    - low-quality sliding window trim from front to tail
    - poly-G tail removal
    - remove reads < min `--length`
    - trim reads down to max `--length`
    """
    workflow = Workflow("qc", "qc.smk", output, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake, no_temp)
    workflow.notebook_files = ["qc_bx_stats.ipynb", "fastp_qc.ipynb"]
    workflow.conda = ["qc"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc = not unlinked, quiet = quiet)
    trim_adapters = False if not trim_adapters else trim_adapters

    workflow.notebooks["skip"] = skip_reports
    workflow.linkedreads["type"] = fastq.lr_type
    workflow.input(fastq.files)
    workflow.param(trim_adapters, "trim-adapters")
    workflow.param(deduplicate, "deduplicate")
    workflow.param(length[0], "min-len")
    workflow.param(length[1], "max-len")
    if extra_params:
        workflow.param(extra_params, "extra")

    treatment = []
    if deduplicate:
        treatment.append("deduplicate")
    if trim_adapters:
        treatment.append("trim adapters")

    workflow.info = {
        "Samples": fastq.count,
        "Linked-Read Type": fastq.lr_type,
        **({"Treatment" : ", ".join(treatment)} if treatment else {"Treatment": "None"}),
        "Output Folder": os.path.relpath(output) + "/",
    }
    
    workflow.initialize(setup)
