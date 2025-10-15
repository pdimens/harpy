"""Harpy demultiplex workflows"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTQfile, DemuxSchema
from harpy.common.cli_types_generic import SnakemakeParams
from harpy.common.printing import workflow_info
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow

@click.group(options_metavar='', context_settings={"help_option_names" : []})
def demultiplex():
    """
    Demultiplex haplotagged FASTQ files

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
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only', panel = "Workflow Options",  is_flag = True, hidden = True, default = False,  help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('schema', required = True, type=DemuxSchema())
@click.argument('R12_FQ', required=True, type=FASTQfile(dir_ok= False), nargs=2)
@click.argument('I12_FQ', required=True, type=FASTQfile(dir_ok= False), nargs=2)
def meier2021(r12_fq, i12_fq, output_dir, schema, qx_rx, keep_unknown_samples, keep_unknown_barcodes, threads, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Demultiplex FASTQ files haplotagged with the Meier _et al._ 2021 protocol

    Use the R1, R2, I2, and I2 FASTQ files provided by the sequencing facility as inputs after the options and schema (4 files, in that exact order). 
    The `SCHEMA` must have **no header** (i.e. no column names) and be in the format of `sample`\\<TAB\\>`barcode`,
    where `barcode` is the barcode segment associated with the sample ID (.e.g. `C01`, `C02`, etc.). Use `--qx-rx` to add the 
    `QX:Z` (barcode PHRED scores) and `RX:Z` (nucleotide barcode) tags in the sequence headers. These tags aren't used by any
    subsequent analyses, but may be useful for your own diagnostics. 
    """
    workflow = Workflow("demultiplex_meier2021", "demultiplex_meier2021.smk", output_dir, quiet) 
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.conda = ["demultiplex", "qc"]

    workflow.inputs = {
        "demultiplex_schema" : schema,
        "R1": r12_fq[0][0],
        "R2": r12_fq[1][0],
        "I1": i12_fq[0][0],
        "I2": i12_fq[1][0]
    }
    workflow.config = {
        "workflow" : workflow.name,
        "retain" : {
            "qx_rx" : qx_rx,
            "barcodes" : keep_unknown_barcodes,
            "samples" : keep_unknown_samples,
        },
        "reports" : {
            "skip": skip_reports
        }
    }
    
    workflow.start_text = workflow_info(
        ("Barcode Design:", "Meier [italic]et al.[/] 2021"),
        ("Demultiplex Schema:", os.path.basename(schema)),
        ("Include QX/RX tags", "Yes" if qx_rx else "No"),
        ("Output Folder:", os.path.relpath(output_dir) + "/")
    )

    workflow.initialize(setup_only)

demultiplex.add_command(meier2021)
