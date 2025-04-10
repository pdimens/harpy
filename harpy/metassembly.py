"""Perform a linked-read aware metassembly"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, instantiate_dir, setup_snakemake, write_workflow_config
from ._cli_types_generic import convert_to_int, HPCProfile, KParam, SnakemakeParams
from ._cli_types_params import SpadesParams
from ._misc import filepath
from ._printing import workflow_info
from ._validations import validate_fastq_bx

docstring = {
    "harpy metassembly": [
        {
            "name": "Metassembly Parameters",
            "options": ["--bx-tag", "--extra-params", "--ignore-bx","--kmer-length", "--max-memory", "--metassembly", "--organism-type"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/metassembly")
# SPADES
@click.option('-b', '--bx-tag', type = click.Choice(['BX', 'BC'], case_sensitive=False), default = "BX", show_default=True, help = "The header tag with the barcode (`BX` or `BC`)")
@click.option('-k', '--kmer-length', type = KParam(), show_default = True, default = "auto", help = 'K values to use for assembly (`odd` and `<128`)')
@click.option('-r', '--max-memory',  type = click.IntRange(min = 1000), show_default = True, default = 10000, help = 'Maximum memory for spades to use, in megabytes')
@click.option('--ignore-bx', is_flag = True, show_default = True, default = False, help = 'Ignore linked-read info for initial spades assembly')
@click.option('-x', '--extra-params', type = SpadesParams(), help = 'Additional spades parameters, in quotes')
# Common Workflow
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Metassembly", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1,999, clamp = True), help = 'Number of threads to use')
@click.option('-u', '--organism-type', type = click.Choice(['prokaryote', 'eukaryote', 'fungus'], case_sensitive=False), default = "eukaryote", show_default=True, help = "Organism type for assembly report [`eukaryote`,`prokaryote`,`fungus`]")
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('fastq_r1', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=1)
@click.argument('fastq_r2', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=1)
def metassembly(fastq_r1, fastq_r2, bx_tag, kmer_length, max_memory, ignore_bx, output_dir, extra_params, container, threads, snakemake, quiet, hpc, organism_type, setup_only, skip_reports):
    """
    Create a metassembly from linked-reads

    The linked-read barcodes must be in `BX:Z` or `BC:Z` FASTQ header tags. If provided, values for `-k` must be
    separated by commas and without spaces (e.g. `-k 15,23,51`). It is strongly recommended to first deconvolve
    the input FASTQ files with `harpy deconvolve`.
    """
    workflow = "metassembly"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    validate_fastq_bx([fastq_r1, fastq_r2], threads, quiet)

    ## setup workflow ##
    command = setup_snakemake(
        workflow,
        "conda" if not container else "conda apptainer",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, f"metassembly.smk")

    conda_envs = ["align", "assembly", "metassembly", "qc", "spades"]
    configs = {
        "workflow" : workflow,
        "snakemake_log" : sm_log,
        "barcode_tag" : bx_tag.upper(),
        "spades" : {
            'ignore_barcodes' : ignore_bx,
            "k" : 'auto' if kmer_length == "auto" else ",".join(map(str,kmer_length)),
            "max_memory" : max_memory,
            **({'extra' : extra_params} if extra_params else {})
        },
        "snakemake_command" : command.rstrip(),
        "conda_environments" : conda_envs,
        "reports" : {
            "skip": skip_reports,
            "organism_type": organism_type
        },
        "inputs": {
            "fastq_r1" : fastq_r1,
            "fastq_r2" : fastq_r2
        }
    }

    write_workflow_config(configs, output_dir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Barcode Tag: ", bx_tag.upper()),
        ("Kmer Length: ", "auto") if kmer_length == "auto" else ("Kmer Length: ", ",".join(map(str,kmer_length))),
        ("Output Folder:", f"{output_dir}/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, workflow, start_text, output_dir, sm_log, quiet, f"workflow/metassembly.summary")