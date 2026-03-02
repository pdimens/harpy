"""Harpy demultiplex workflows"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTQfile, DemuxSchema
from harpy.common.cli_params import SnakemakeParams
from harpy.common.file_ops import fetch_template
from harpy.common.system_ops import container_ok
from harpy.validation.fastq import FASTQ
from harpy.common.workflow import Workflow

#TODO update docs link to dedicated pages
@click.group(options_metavar='')
@click.help_option('--help', hidden = True)
def preprocess():
    """
    Remove inline barcodes from raw FASTQs

    The provided methods are specific to Haplotagging-style linked-read sequencing.
    Check that you are using the correct haplotagging method/technology, since the different
    barcoding approaches have very different demultiplexing strategies.

    **Haplotagging Technologies**
    - `meier2021`: the original haplotagging barcode strategy
      - Meier _et al._ (2021) doi: 10.1073/pnas.2015005118
    - `gih`: updated/modified protocol developed by the Cornell GIH
      - Iqbal _el al._ (in prep)
    """

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.option('-u', '--keep-unknown-samples', panel = "Parameters",  is_flag = True, default = False, help = 'Keep a separate file of reads with recognized barcodes but don\'t match any sample in the schema')
@click.option('-b', '--keep-unknown-barcodes', panel = "Parameters",  is_flag = True, default = False, help = 'Keep a separate file of reads with unrecognized barcodes')
@click.option('-q', '--qx-rx', panel = "Parameters", is_flag = True, default = False, help = 'Include the `QX:Z` and `RX:Z` tags in the read header')
@click.option('-@', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(2,999, clamp = True), help = 'Number of threads to use')
@click.option('-O', '--output', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Preprocess", show_default=True,  help = 'Output directory name')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('-C', '--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('-N', '--setup', panel = "Workflow Options",  is_flag = True, hidden = True, default = False,  help = 'Setup the workflow and exit')
@click.option('-H', '--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('-Q', '--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('-R', '--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('-S', '--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('schema', required = True, type=DemuxSchema())
@click.argument('R12_FQ', required=True, type=FASTQfile(dir_ok= False), nargs=2)
@click.argument('I12_FQ', required=True, type=FASTQfile(dir_ok= False), nargs=2)
@click.help_option('--help', hidden = True)
def meier2021(r12_fq, i12_fq, output, schema, qx_rx, keep_unknown_samples, keep_unknown_barcodes, threads, snakemake, skip_reports, quiet, hpc, clean, container, setup):
    """
    Preprocess FASTQ files haplotagged with the Meier _et al._ 2021 protocol

    Use the R1, R2, I2, and I2 FASTQ files provided by the sequencing facility as inputs after the options and schema (4 files, in that exact order). 
    The `SCHEMA` must have **no header** (i.e. no column names) and be in the format of `sample`\\<TAB\\>`barcode`,
    where `barcode` is the barcode segment associated with the sample ID (.e.g. `C01`, `C02`, etc.). Use `--qx-rx` to add the 
    `QX:Z` (barcode PHRED scores) and `RX:Z` (nucleotide barcode) tags in the sequence headers. These tags aren't used by any
    subsequent analyses, but may be useful for your own diagnostics. 
    """
    workflow = Workflow("preprocess_meier2021", "preprocess_meier2021.smk", output, container, clean, quiet, no_validation=True) 
    workflow.setup_snakemake(threads, hpc, snakemake, no_temp)
    workflow.conda = ["qc"]
    
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
    
    workflow.info = {
        "Barcode Design": "Meier [italic]et al.[/] 2021",
        "Demultiplex Schema": os.path.basename(schema),
        "Include QX/RX tags" : "Yes" if qx_rx else "No",
        "Output Folder" : os.path.relpath(output) + "/"
    }

    workflow.initialize(setup)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.option('-m', '--me-seq', panel = "Parameters", type = str, default = "AGATGTGTATAAGAGACAG", show_default=True, help = "ME sequence to look for")
@click.option('-o', '--me-overlap', panel = "Parameters", default = 19, show_default = True, type = click.IntRange(0,300, clamp = True), help = 'ME sequence overlap')
@click.option('-@', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(2,999, clamp = True), help = 'Number of threads to use')
@click.option('-O', '--output', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Preprocess", show_default=True,  help = 'Output directory name')
@click.option('-T', '--no-temp', hidden = True, panel = "Workflow Options", is_flag = True, default = False, help = 'Don\'t delete temporary files')
@click.option('-C', '--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('-N', '--setup', panel = "Workflow Options",  is_flag = True, hidden = True, default = False,  help = 'Setup the workflow and exit')
@click.option('-H', '--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('-Q', '--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('-R', '--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('-S', '--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
@click.help_option('--help', hidden = True)
def gih(inputs, output, me_seq, me_overlap, threads, snakemake, skip_reports, quiet, hpc, clean, container, setup, no_temp):
    """
    Preprocess FASTQ files haplotagged with the GIH protocol

    
    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/poccidentalis*.fq`), or both.
    The resulting FASTQ file pairs will have inline barcodes removed and added to the sequence headers,
    but will still need QC to remove adapters and low-quality reads/regions. 
    """
    workflow = Workflow("preprocess_gih", "preprocess_gih.smk", output, container, clean, quiet, no_validation=True) 
    workflow.setup_snakemake(threads, hpc, snakemake, no_temp)
    workflow.conda = ["qc", "preprocess"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc = False, quiet= True)

    workflow.notebooks["skip"] = skip_reports
    fetch_template("pheniqs.config.json", os.path.join(output, "workflow", "pheniqs.config.json"))
    workflow.input(fastq.files)
    workflow.param(me_seq, "ME-sequence")
    workflow.param(me_overlap, "ME-overlap")
    
    workflow.info = {
        "Barcode Design": "Iqbal [italic]et al.[/] (in prep)",
        "Output Folder": os.path.relpath(output) + "/"
    }

    workflow.initialize(setup)

preprocess.add_command(meier2021)
preprocess.add_command(gih)
