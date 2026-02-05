"""Harpy demultiplex workflows"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTQfile, DemuxSchema
from harpy.common.cli_types_generic import SnakemakeParams
from harpy.common.printing import workflow_info
from harpy.common.system_ops import container_ok
from harpy.validation.fastq import FASTQ
from harpy.common.workflow import Workflow

@click.group(options_metavar='')
@click.help_option('--help', hidden = True)
def preprocess():
    """
    preprocess haplotagging FASTQ files

    Check that you are using the correct haplotagging method/technology, since the different
    barcoding approaches have very different demultiplexing strategies.

    **Haplotagging Technologies**
    - `meier2021`: the original haplotagging barcode strategy
      - Meier _et al._ (2021) doi: 10.1073/pnas.2015005118
    """

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/demultiplex/")
@click.option('-u', '--keep-unknown-samples', panel = "Parameters",  is_flag = True, default = False, help = 'Keep a separate file of reads with recognized barcodes but don\'t match any sample in the schema')
@click.option('-b', '--keep-unknown-barcodes', panel = "Parameters",  is_flag = True, default = False, help = 'Keep a separate file of reads with unrecognized barcodes')
@click.option('-q', '--qx-rx', panel = "Parameters", is_flag = True, default = False, help = 'Include the `QX:Z` and `RX:Z` tags in the read header')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(2,999, clamp = True), help = 'Number of threads to use')
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Demultiplex", show_default=True,  help = 'Output directory name')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup', panel = "Workflow Options",  is_flag = True, hidden = True, default = False,  help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('schema', required = True, type=DemuxSchema())
@click.argument('R12_FQ', required=True, type=FASTQfile(dir_ok= False), nargs=2)
@click.argument('I12_FQ', required=True, type=FASTQfile(dir_ok= False), nargs=2)
@click.help_option('--help', panel = "Workflow Options", hidden = True)
def meier2021(r12_fq, i12_fq, output_dir, schema, qx_rx, keep_unknown_samples, keep_unknown_barcodes, threads, snakemake, skip_reports, quiet, hpc, clean, container, setup):
    """
    Demultiplex FASTQ files haplotagged with the Meier _et al._ 2021 protocol

    Use the R1, R2, I2, and I2 FASTQ files provided by the sequencing facility as inputs after the options and schema (4 files, in that exact order). 
    The `SCHEMA` must have **no header** (i.e. no column names) and be in the format of `sample`\\<TAB\\>`barcode`,
    where `barcode` is the barcode segment associated with the sample ID (.e.g. `C01`, `C02`, etc.). Use `--qx-rx` to add the 
    `QX:Z` (barcode PHRED scores) and `RX:Z` (nucleotide barcode) tags in the sequence headers. These tags aren't used by any
    subsequent analyses, but may be useful for your own diagnostics. 
    """
    workflow = Workflow("demultiplex_meier2021", "demultiplex_meier2021.smk", output_dir, container, clean, quiet) 
    workflow.setup_snakemake(threads, hpc, snakemake)
    workflow.conda = ["demultiplex", "qc"]
    
    workflow.inputs = {
        "schema" : schema,
        "R1": r12_fq[0][0],
        "R2": r12_fq[1][0],
        "I1": i12_fq[0][0],
        "I2": i12_fq[1][0]
    }
    workflow.param("qx-rx", qx_rx)
    workflow.param("barcodes", keep_unknown_barcodes)
    workflow.param("samples", keep_unknown_samples)
    workflow.notebooks["skip"] = skip_reports
    
    workflow.start_text = workflow_info(
        ("Barcode Design:", "Meier [italic]et al.[/] 2021"),
        ("Demultiplex Schema:", os.path.basename(schema)),
        ("Include QX/RX tags", "Yes" if qx_rx else "No"),
        ("Output Folder:", os.path.relpath(output_dir) + "/")
    )

    workflow.initialize(setup)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/demultiplex/")
@click.option('-l', '--spacer-length', panel = "Parameters",  type = click.IntRange(min = 10), default = 77, help = 'Length of spacers between barcodes')
@click.option('-m', '--min-length', panel = "Parameters", type = click.IntRange(min = 5),  default = 50, help = 'Minimum insert length (bp) of reads to retain')
@click.option('-q', '--min-quality', panel = "Parameters", type = click.IntRange(min = 0), default = 20, help = 'Minimum average read quality to retain')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(2,999, clamp = True), help = 'Number of threads to use')
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Demultiplex", show_default=True,  help = 'Output directory name')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup', panel = "Workflow Options",  is_flag = True, hidden = True, default = False,  help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('barcodes', required = True, type=DemuxSchema())
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
@click.help_option('--help', panel = "Workflow Options", hidden = True)
def gih(inputs, output_dir, barcodes, spacer_length, min_length, min_quality, threads, snakemake, skip_reports, quiet, hpc, clean, container, setup):
    """
    Demultiplex FASTQ files haplotagged with the Genomics Innovation Hub protocol

    
    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/poccidentalis*.fq`), or both.
    The `BARCODES` file must have **no header** (i.e. no column name). 
    """
    workflow = Workflow("demultiplex_gih", "demultiplex_gih.smk", output_dir, container, clean, quiet) 
    workflow.setup_snakemake(threads, hpc, snakemake)
    workflow.conda = ["demultiplex", "qc"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc = False, quiet= True)
    with open(barcodes, "r") as f:
        bc_seg_len = len(f.readline())

    bc_len = (3 * bc_seg_len) + (2 * spacer_length)
    bc_len_text = f"{bc_len} (3×barcode + 2×spacer)"

    workflow.notebooks["skip"] = skip_reports
    workflow.input(fastq.files)
    workflow.param(bc_len, "barcode_length")
    workflow.param(min_length, "minimum_length")
    workflow.param(min_quality, "minimum_quality")
    
    workflow.start_text = workflow_info(
        ("Barcode Design:", "Iqbal [italic]et al.[/] (in prep)"),
        ("Total Barcode Length:", bc_len_text),
        ("Min. insert length:", min_length),
        ("Min. read quality:", min_quality),
        ("Output Folder:", os.path.relpath(output_dir) + "/")
    )

    workflow.initialize(setup)

preprocess.add_command(meier2021)
preprocess.add_command(gih)
