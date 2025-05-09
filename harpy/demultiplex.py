"""Harpy demultiplex workflows"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._cli_types_generic import convert_to_int, HPCProfile, SnakemakeParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, instantiate_dir, setup_snakemake, write_workflow_config
from ._printing import workflow_info
from ._validations import validate_demuxschema

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def demultiplex():
    """
    Demultiplex haplotagged FASTQ files

    Check that you are using the correct haplotagging method/technology, since the different
    barcoding approaches have very different demultiplexing strategies.

    **Haplotagging Technologies**
    - `meier2021`: the original haplotagging barcode strategy
      - Meier _et al._ (2021) doi: 10.1073/pnas.2015005118
    """

docstring = {
    "harpy demultiplex meier2021": [
        {
            "name": "Parameters",
            "options": ["--keep-unknown-barcodes", "--keep-unknown-samples", "--qx-rx","--schema"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(deprecated = True)
def gen1():
    print("This workflow has been renamed \"meier2021\"-- please use that instead. This warning will be removed in the next minor Harpy version and will only return an error.")
    sys.exit(1)

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/demultiplex/")
@click.option('-u', '--keep-unknown-samples',  is_flag = True, default = False, help = 'Keep a separate file of reads with recognized barcodes but don\'t match any sample in the schema')
@click.option('-b', '--keep-unknown-barcodes',  is_flag = True, default = False, help = 'Keep a separate file of reads with unrecognized barcodes')
@click.option('-q', '--qx-rx', is_flag = True, default = False, help = 'Include the `QX:Z` and `RX:Z` tags in the read header')
@click.option('-s', '--schema', required = True, type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True), help = 'File of `sample`\\<TAB\\>`barcode`')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(2,999, clamp = True), help = 'Number of threads to use')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Demultiplex", show_default=True,  help = 'Output directory name')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False,  help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('R1_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True))
@click.argument('R2_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True))
@click.argument('I1_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True))
@click.argument('I2_FQ', required=True, type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True))
def meier2021(r1_fq, r2_fq, i1_fq, i2_fq, output_dir, schema, qx_rx, keep_unknown_samples, keep_unknown_barcodes, threads, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Demultiplex FASTQ files haplotagged with the Meier _et al._ 2021 protocol

    Use the R1, R2, I2, and I2 FASTQ files provided by the sequencing facility as inputs (in that exact order) provided after the options. 
    The `--schema` must be **tab** (or space) delimited, have **no header** (i.e. no column names), and be in the format of `sample`\\<TAB\\>`barcode`,
    where `barcode` is the barcode segment associated with the sample ID (.e.g. `C01`, `C02`, etc.). Use `--qx-rx` to add the 
    `QX:Z` (barcode PHRED scores) and `RX:Z` (nucleotide barcode) tags in the sequence headers. These tags aren't used by any
    subsequent analyses, but may be useful for your own diagnostics. 
    """
    workflow = "demultiplex_meier2021"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    validate_demuxschema(schema, return_len = False)

    ## setup workflow ##
    command,command_rel = setup_snakemake(
        workflow,
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "demultiplex_meier2021.smk")

    conda_envs = ["demultiplex", "qc"]
    configs = {
        "workflow" : workflow,
        "retain" : {
            "qx_rx" : qx_rx,
            "barcodes" : keep_unknown_barcodes,
            "samples" : keep_unknown_samples,
        },
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports
        },
        "inputs" : {
            "demultiplex_schema" : schema,
            "R1": r1_fq,
            "R2": r2_fq,
            "I1": i1_fq,
            "I2": i2_fq
        }
    }
    
    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)
    
    start_text = workflow_info(
        ("Barcode Design:", "Meier _et al._ 2021"),
        ("Demultiplex Schema:", os.path.basename(schema)),
        ("Include QX/RX tags", "Yes" if qx_rx else "No"),
        ("Output Folder:", os.path.basename(output_dir) + "/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, "workflow/demux.meier2021.summary")

demultiplex.add_command(meier2021)
demultiplex.add_command(gen1)
