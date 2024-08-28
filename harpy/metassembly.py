"""Perform a linked-read aware metassembly"""

import os
import sys
from rich import box
from rich.table import Table
import rich_click as click
from ._conda import generate_conda_deps
from ._launch import launch_snakemake
from ._misc import fetch_rule, snakemake_log, IntPair

docstring = {
    "harpy metassembly": [
        {
            "name": "Parameters",
            "options": ["--clusters", "--contig-cov"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--skipreports", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/qc")
@click.option('-n', '--clusters', default = 35, show_default = True, type = int, help = 'Number of clusters')
@click.option('-c', '--contig-cov', default = "10,30", show_default = True, type = IntPair(), help = "Coverage for low abundance contigs")
@click.option('-x', '--extra-params', type = str, help = 'Additional STITCH parameters, in quotes')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Metassembly", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 1, max_open = True), help = 'Number of threads to use')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake',  type = str, help = 'Additional Snakemake parameters, in quotes')
@click.argument('fastq', required=True, type=click.Path(exists=True, readable=True), nargs=1)
@click.argument('fastq2', required=False, type=click.Path(exists=True, readable=True), nargs=1)
def metassembly(fastq, fastq2, clusters, contig_cov, output_dir, extra_params, threads, snakemake, skipreports, quiet, hpc, conda, setup_only):
    """
    Do a metassembly
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock  --software-deployment-method {sdm} --conda-prefix .snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/metassembly.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake is not None:
        command += snakemake

    os.makedirs(workflowdir, exist_ok=True)
    fetch_rule(workflowdir, "metassembly.smk")
    #fetch_report(workflowdir, "bx_count.Rmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "metassembly")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: metassembly\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"clusters: {clusters}\n")
        config.write(f"contig_coverage: {contig_cov[0]},{contig_cov[1]}\n")
        if extra_params:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skip_reports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  fastq: {fastq}\n")
        if fastq2:
            config.write(f"  fastq2: {fastq2}\n")

    generate_conda_deps()
    if setup_only:
        sys.exit(0)

    n_fq = 2 if fastq2 else 1
    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column("value", justify="left")
    start_text.add_row("FASTQ Inputs: ", f"{n_fq}")
    start_text.add_row("Clusters: ", f"{clusters}")
    start_text.add_row("Contig Cov. Thresh: ", f"{contig_cov[0]},{contig_cov[1]}")
    start_text.add_row("Output Folder:", f"{output_dir}/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", ""))
    launch_snakemake(command, "metassembly", start_text, output_dir, sm_log, quiet)
