"""Perform a linked-read aware metassembly"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._cli_types_generic import convert_to_int, KParam, HPCProfile, SnakemakeParams
from ._cli_types_params import SpadesParams, ArcsParams
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_rule, snakemake_log, write_snakemake_config, write_workflow_config
from ._printing import workflow_info
from ._validations import validate_fastq_bx

docstring = {
    "harpy assembly": [
        {
            "name": "Assembly Parameters",
            "options": ["--bx-tag", "--extra-params", "--kmer-length", "--max-memory", "--organism-type"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Scaffolding Parameters",
            "options": ["--arcs-extra", "--contig-length", "--links", "--min-aligned", "--min-quality", "--mismatch", "--molecule-distance", "--molecule-length", "--seq-identity", "--span"],            
            "panel_styles": {"border_style": "green"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/assembly")
# SPADES
@click.option('-b', '--bx-tag', type = click.Choice(['BX', 'BC'], case_sensitive=False), default = "BX", show_default=True, help = "The header tag with the barcode [`BX`,`BC`]")
@click.option('-k', '--kmer-length', type = KParam(), show_default = True, default = "auto", help = 'K values to use for assembly (`odd` and `<128`)')
@click.option('-r', '--max-memory',  type = click.IntRange(min = 1000), show_default = True, default = 10000, help = 'Maximum memory for spades to use, in megabytes')
@click.option('-x', '--extra-params', type = SpadesParams(), help = 'Additional spades parameters, in quotes')
# TIGMINT/ARCS/LINKS
@click.option('-y', '--arcs-extra', type = ArcsParams(), help = 'Additional ARCS parameters, in quotes (`option=arg` format)')
@click.option("-c","--contig-length", type = int, default = 500, show_default = True, help = "Minimum contig length")
@click.option("-n", "--links", type = int, default = 5, show_default = True, help = "Minimum number of links to compute scaffold")
@click.option("-a", "--min-aligned", type = int, default = 5, show_default = True, help = "Minimum aligned read pairs per barcode")
@click.option("-q", "--min-quality", type = click.IntRange(0,40, clamp = True), default = 0, show_default = True, help = "Minimum mapping quality")
@click.option("-m", "--mismatch", type = int, default = 5, show_default = True, help = "Maximum number of mismatches")
@click.option("-d", "--molecule-distance", type = int, default = 50000, show_default = True, help = "Distance cutoff to split molecules (bp)")
@click.option("-l", "--molecule-length", type = int, default = 2000, show_default = True, help = "Minimum molecule length (bp)")
@click.option("-i", "--seq-identity", type = click.IntRange(0,100, clamp = True), default = 98, show_default = True, help = "Minimum sequence identity") 
@click.option("-s", "--span", type = int, default = 20, show_default = True, help = "Minimum number of spanning molecules to be considered assembled")
# Other Options
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Assembly", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1, 999, clamp = True), help = 'Number of threads to use')
@click.option('-u', '--organism-type', type = click.Choice(['prokaryote', 'eukaryote', 'fungus'], case_sensitive=False), default = "eukaryote", show_default=True, help = "Organism type for assembly report [`eukaryote`,`prokaryote`,`fungus`]")
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('fastq_r1', required=True, type=click.Path(exists=True, readable=True), nargs=1)
@click.argument('fastq_r2', required=True, type=click.Path(exists=True, readable=True), nargs=1)
def assembly(fastq_r1, fastq_r2, bx_tag, kmer_length, max_memory, output_dir, extra_params,arcs_extra,contig_length,links,min_quality,min_aligned,mismatch,molecule_distance,molecule_length,seq_identity,span, organism_type, container, threads, snakemake, quiet, hpc, setup_only, skip_reports):
    """
    Create an assembly from linked-reads

    The linked-read barcodes must be in `BX:Z` or `BC:Z` FASTQ header tags. If provided, values for `-k` must be
    separated by commas and without spaces (e.g. `-k 15,23,51`). It is strongly recommended to first deconvolve
    the input FASTQ files with `harpy deconvolve`.
    """
    ## checks and validations ##
    validate_fastq_bx([fastq_r1, fastq_r2], threads, quiet)

    ## setup workflow #
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    write_snakemake_config("conda" if not container else "conda apptainer", output_dir)
    command = f"snakemake --cores {threads} --snakefile {workflowdir}/assembly.smk"
    command += f" --configfile {workflowdir}/workflow.yaml --profile {workflowdir}"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    fetch_rule(workflowdir, f"assembly.smk")

    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "assembly")
    conda_envs = ["assembly","qc"]
    configs = {
        "workflow" : "assembly",
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
        "workflow_call" : command.rstrip(),
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

    write_workflow_config(configs, workflowdir)
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Barcode Tag: ", bx_tag.upper()),
        ("Kmer Length: ", "auto") if kmer_length == "auto" else ("Kmer Length: ", ",".join(map(str,kmer_length))),
        ("Output Folder:", f"{output_dir}/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "assembly", start_text, output_dir, sm_log, quiet, f"workflow/assembly.summary")