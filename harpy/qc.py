"""Harpy sequence adapter trimming and quality control"""

import os
import sys
from rich import box
from rich.table import Table
import rich_click as click
from ._conda import create_conda_recipes
from ._launch import launch_snakemake
from ._misc import fetch_report, fetch_rule, snakemake_log, IntQuartet
from ._parsers import parse_fastq_inputs

docstring = {
    "harpy qc": [
        {
            "name": "Parameters",
            "options": ["--deconvolve", "--deconvolve-params", "--deduplicate", "--extra-params", "--min-length", "--max-length", "--trim-adapters"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "See the documentation for more information: https://pdimens.github.io/harpy/workflows/qc")
@click.option('-c', '--deconvolve', is_flag = True, default = False, help = 'Resolve barcode clashes between reads from different molecules.')
@click.option('-p', '--deconvolve-params', type = IntQuartet(), show_default = True, default = "21,40,3,0", help = 'Accepts the QuickDeconvolution parameters for k,w,d,a, in that order')
@click.option('-d', '--deduplicate', is_flag = True, default = False, help = 'Identify and remove PCR duplicates')
@click.option('-x', '--extra-params', type = str, help = 'Additional Fastp parameters, in quotes')
@click.option('-m', '--max-length', default = 150, show_default = True, type=int, help = 'Maximum length to trim sequences down to')
@click.option('-n', '--min-length', default = 30, show_default = True, type=int, help = 'Discard reads shorter than this length')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "QC", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-a', '--trim-adapters', is_flag = True, default = False, help = 'Detect and trim adapters')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skip-reports',  is_flag = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def qc(inputs, output_dir, min_length, max_length, trim_adapters, deduplicate, deconvolve, deconvolve_params, extra_params, threads, snakemake, skip_reports, quiet, hpc, conda, setup_only):
    """
    Remove adapters and quality-control sequences

    Provide the input fastq files and/or directories at the end of the command
    as individual files/folders, using shell wildcards (e.g. `data/acronotus*.fq`), or both.
    
    The input reads will be quality trimmed using:
    - a sliding window from front to tail
    - poly-G tail removal
    - use `-a` to automatically detect and remove adapters
    - use `-d` to find and remove PCR duplicates
    - use `-c` to resolve barcode clashing that may occur by unrelated sequences having the same barcode
      - the parameters for `-p` are described [here](https://github.com/RolandFaure/QuickDeconvolution?tab=readme-ov-file#usage).
      - you can use `harpy deconvolve` to perform this task separately
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock  --software-deployment-method {sdm} --conda-prefix .snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/qc.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake is not None:
        command += snakemake

    fqlist, sample_count = parse_fastq_inputs(inputs)
    os.makedirs(workflowdir, exist_ok=True)
    fetch_rule(workflowdir, "qc.smk")
    fetch_report(workflowdir, "bx_count.Rmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "qc")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: qc\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"trim_adapters: {trim_adapters}\n")
        config.write(f"deduplicate: {deduplicate}\n")
        if deconvolve:
            config.write("deconvolve:\n")
            k,w,d,a = deconvolve_params
            config.write(f"  kmer_length: {k}\n")
            config.write(f"  window_size: {w}\n")
            config.write(f"  density: {d}\n")
            config.write(f"  dropout: {a}\n")
        config.write(f"min_len: {min_length}\n")
        config.write(f"max_len: {max_length}\n")
        if extra_params:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skip_reports: {skip_reports}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        for i in fqlist:
            config.write(f"  - {i}\n")

    create_conda_recipes()
    if setup_only:
        sys.exit(0)

    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column("value", justify="left")
    start_text.add_row("Samples:", f"{sample_count}")
    start_text.add_row("Trim Adapters:", "yes" if trim_adapters else "no")
    start_text.add_row("Deduplicate:", "yes" if deduplicate else "no")
    start_text.add_row("Deconvolve:", "yes" if deconvolve else "no")
    start_text.add_row("Output Folder:", f"{output_dir}/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    launch_snakemake(command, "qc", start_text, output_dir, sm_log, quiet)
