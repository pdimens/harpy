"""Perform a linked-read aware metassembly"""

import os
import sys
import yaml
from rich import box
from rich.table import Table
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, snakemake_log, KParam
from ._validations import validate_fastq_bx

docstring = {
    "harpy assembly": [
        {
            "name": "Assembly Parameters",
            "options": ["--bx-tag", "--extra-params", "--kmer-length", "--max-memory"],
        },
        {
            "name": "Scaffolding Parameters (ignored for metassembly)",
            "options": ["--arcs-extra", "--contig-length", "--links", "--min-aligned", "--min-quality", "--mismatch", "--molecule-distance", "--molecule-length", "--seq-identity", "--span"],            
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/assembly")
# SPADES
@click.option('-b', '--bx-tag', type = click.Choice(['BX', 'BC'], case_sensitive=False), default = "BX", show_default=True, help = "The header tag with the barcode (`BX` or `BC`)")
@click.option('-k', '--kmer-length', type = KParam(), show_default = True, default = "auto", help = 'K values to use for assembly (`odd` and `<128`)')
@click.option('-r', '--max-memory',  type = click.IntRange(min = 1000, max_open = True), show_default = True, default = 10000, help = 'Maximum memory for spades to use, in megabytes')
@click.option('-x', '--extra-params', type = str, help = 'Additional spades parameters, in quotes')
# TIGMINT/ARCS/LINKS
@click.option('-y', '--arcs-extra', type = str, help = 'Additional ARCS parameters, in quotes')
@click.option("-c","--contig-length", type = int, default = 500, show_default = True, help = "Minimum contig length")
@click.option("-n", "--links", type = int, default = 5, show_default = True, help = "Minimum number of links to compute scaffold")
@click.option("-a", "--min-aligned", type = int, default = 5, show_default = True, help = "Minimum aligned read pairs per barcode")
@click.option("-q", "--min-quality", type = click.IntRange(0,40), default = 0, show_default = True, help = "Minimum mapping quality")
@click.option("-m", "--mismatch", type = int, default = 5, show_default = True, help = "Maximum number of mismatches")
@click.option("-d", "--molecule-distance", type = int, default = 50000, show_default = True, help = "Distance cutoff to split molecules (bp)")
@click.option("-l", "--molecule-length", type = int, default = 2000, show_default = True, help = "Minimum molecule length (bp)")
@click.option("-i", "--seq-identity", type = click.IntRange(0,100), default = 98, show_default = True, help = "Minimum sequence identity") 
@click.option("-s", "--span", type = int, default = 20, show_default = True, help = "Minimum number of spanning molecules to be considered assembled")
# Common Workflow
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Assembly", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake',  type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('fastq_r1', required=True, type=click.Path(exists=True, readable=True), nargs=1)
@click.argument('fastq_r2', required=True, type=click.Path(exists=True, readable=True), nargs=1)
def assembly(fastq_r1, fastq_r2, bx_tag, kmer_length, max_memory, output_dir, extra_params,arcs_extra,contig_length,links,min_quality,min_aligned,mismatch,molecule_distance,molecule_length,seq_identity,span,conda, threads, snakemake, skip_reports, quiet, hpc, setup_only):
    """
    Create an assembly from linked-reads

    The linked-read barcodes must be in `BX:Z` or `BC:Z` FASTQ header tags. If provided, values for `-k` must be
    separated by commas and without spaces (e.g. `-k 15,23,51`). It is strongly recommended to first deconvolve
    the input FASTQ files with `harpy deconvolve`.
    """
    output_dir = output_dir.rstrip("/")
    asm = "assembly"
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock  --software-deployment-method {sdm} --conda-prefix .snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/{asm}.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake:
        command += snakemake

    validate_fastq_bx([fastq_r1, fastq_r2], threads, quiet)
    os.makedirs(workflowdir, exist_ok=True)
    fetch_rule(workflowdir, f"{asm}.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, asm)
    configs = {
        "workflow" : asm,
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "barcode_tag" : bx_tag.upper(),
        "spades" : {
            "k" : 'auto' if kmer_length == "auto" else ",".join(map(str,kmer_length)),
            "max_memory" : max_memory,
            **({'extra' : extra_params} if extra_params else {})
        },
        "tigmint" : {
            "minimum_mapping_quality" : min_quality,
            "mismatch" : mismatch,
            "molecule_distance" : molecule_distance,
            "molecule_length" : molecule_length,
            "span" : span
        },
        "arcs" : {
            "minimum_aligned_reads" : min_aligned,
            "minimum_contig_length" : contig_length,
            "minimum_sequence_identity" : seq_identity,
            **({'extra' : arcs_extra} if arcs_extra else {})
        },
        "links" : {
            "minimum_links" : links
        },
        "skip_reports" : skip_reports,
        "workflow_call" : command.rstrip(),
        "inputs": {
            "fastq_r1" : fastq_r1,
            "fastq_r2" : fastq_r2
        }
    }
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))
    
    create_conda_recipes()
    if setup_only:
        sys.exit(0)

    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column("value", justify="left")
    start_text.add_row("Barcode Tag: ", bx_tag.upper())
    if kmer_length == "auto":
        start_text.add_row("Kmer Length: ", "auto")
    else:
        start_text.add_row("Kmer Length: ", ",".join(map(str,kmer_length)))
    start_text.add_row("Output Folder:", f"{output_dir}/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    launch_snakemake(command, asm, start_text, output_dir, sm_log, quiet, f"workflow/{asm}.summary")