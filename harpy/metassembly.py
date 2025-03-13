"""Perform a linked-read aware metassembly"""

import os
import sys
import yaml
import shutil
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake, SNAKEMAKE_CMD
from ._misc import fetch_rule, snakemake_log
from ._cli_types_generic import convert_to_int, HPCProfile, KParam, SnakemakeParams
from ._cli_types_params import SpadesParams
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
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Metassembly", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1,999, clamp = True), help = 'Number of threads to use')
@click.option('-u', '--organism-type', type = click.Choice(['prokaryote', 'eukaryote', 'fungus'], case_sensitive=False), default = "eukaryote", show_default=True, help = "Organism type for assembly report [`eukaryote`,`prokaryote`,`fungus`]")
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('fastq_r1', required=True, type=click.Path(exists=True, readable=True), nargs=1)
@click.argument('fastq_r2', required=True, type=click.Path(exists=True, readable=True), nargs=1)
def metassembly(fastq_r1, fastq_r2, bx_tag, kmer_length, max_memory, ignore_bx, output_dir, extra_params, container, threads, snakemake, quiet, hpc, organism_type, setup_only, skip_reports):
    """
    Set up and optionally run a linked-read metassembly workflow.
    
    This function configures and executes a Snakemake-based metassembly workflow for linked-read data.
    It validates paired FASTQ files with barcode headers (in "BX:Z" or "BC:Z" format), prepares the workflow
    directory, writes the configuration (including SPAdes parameters and optional HPC settings), and creates
    the required conda environments. When the 'setup_only' flag is True, the workflow is configured and the
    process exits without launching Snakemake.
    
    Args:
        fastq_r1: Path to the first FASTQ file.
        fastq_r2: Path to the second FASTQ file.
        bx_tag: Barcode tag identifier (e.g. "BX" or "BC") used to extract linked-read barcodes.
        kmer_length: Kmer lengths for SPAdes assembly as a list of integers or "auto" for automatic selection.
        max_memory: Maximum memory allocation for SPAdes.
        ignore_bx: Flag indicating whether to ignore linked-read barcodes during assembly.
        output_dir: Directory where workflow files, logs, and outputs will be stored.
        extra_params: Additional parameters for SPAdes assembly.
        container: Boolean flag indicating whether to deploy software using a container.
        threads: Number of threads to allocate for the workflow.
        snakemake: Extra command-line arguments to pass to Snakemake.
        quiet: Flag to suppress verbose output.
        hpc: Path to an HPC submission YAML configuration file; if provided, it is copied to the workflow's
             HPC directory.
        organism_type: Specifies the organism type for report configuration.
        setup_only: If True, only sets up the workflow without launching Snakemake.
        skip_reports: If True, generation of reports is skipped.
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = os.path.join(output_dir, 'workflow')
    sdm = "conda" if not container else "conda apptainer"
    command = f'{SNAKEMAKE_CMD} --software-deployment-method {sdm} --cores {threads}'
    command += f" --snakefile {workflowdir}/metassembly.smk"
    command += f" --configfile {workflowdir}/config.yaml"
    if hpc:
        os.makedirs(f"{workflowdir}/hpc", exist_ok=True)
        shutil.copy2(hpc, f"{workflowdir}/hpc/config.yaml")
        command += f" --workflow-profile {workflowdir}/hpc"
    if snakemake:
        command += f" {snakemake}"

    validate_fastq_bx([fastq_r1, fastq_r2], threads, quiet)
    os.makedirs(workflowdir, exist_ok=True)
    fetch_rule(workflowdir, f"metassembly.smk")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "metassembly")
    conda_envs = ["align", "assembly", "metassembly", "qc", "spades"]
    configs = {
        "workflow" : "metassembly",
        "snakemake_log" : sm_log,
        "output_directory" : output_dir,
        "barcode_tag" : bx_tag.upper(),
        "spades" : {
            'ignore_barcodes' : ignore_bx,
            "k" : 'auto' if kmer_length == "auto" else ",".join(map(str,kmer_length)),
            "max_memory" : max_memory,
            **({'extra' : extra_params} if extra_params else {})
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
    with open(os.path.join(workflowdir, 'config.yaml'), "w", encoding="utf-8") as config:
        yaml.dump(configs, config, default_flow_style= False, sort_keys=False, width=float('inf'))
    
    create_conda_recipes(output_dir, conda_envs)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Barcode Tag: ", bx_tag.upper()),
        ("Kmer Length: ", "auto") if kmer_length == "auto" else ("Kmer Length: ", ",".join(map(str,kmer_length))),
        ("Output Folder:", f"{output_dir}/"),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command, "metassembly", start_text, output_dir, sm_log, quiet, f"workflow/metassembly.summary")