"""Harpy imputation workflow"""

import os
import sys
from pathlib import Path
from rich import box
from rich.table import Table
import rich_click as click
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_rule, fetch_report, fetch_script, snakemake_log, launch_snakemake
from .fileparsers import parse_alignment_inputs, biallelic_contigs
from .validations import validate_input_by_ext, vcf_samplematch, check_impute_params, validate_bam_RG

docstring = {
        "harpy impute": [
        {
            "name": "Parameters",
            "options": ["--extra-params",  "--parameters", "--vcf", "--vcf-samples"],
        },
        {
            "name": "Workflow Controls",
            "options": ["--conda", "--hpc", "--output-dir", "--quiet", "--skipreports", "--snakemake", "--threads", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog = "See the documentation for more information: https://pdimens.github.io/harpy/modules/impute/")
@click.option('-x', '--extra-params', type = str, help = 'Additional STITCH parameters, in quotes')
@click.option('-o', '--output-dir', type = click.Path(exists = False), default = "Impute", show_default=True,  help = 'Output directory name')
@click.option('-p', '--parameters', required = True, type=click.Path(exists=True, dir_okay=False), help = 'STITCH parameter file (tab-delimited)')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-v', '--vcf', required = True, type=click.Path(exists=True, dir_okay=False, readable=True), help = 'Path to BCF/VCF file')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False, readable=True), help = 'Directory with HPC submission `config.yaml` file')
@click.option('--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for imputation rather than those found the inputs')
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True), nargs=-1)
def impute(inputs, output_dir, parameters, threads, vcf, vcf_samples, extra_params, snakemake, skipreports, quiet, hpc, conda, setup_only):
    """
    Impute genotypes using variants and alignments
    
    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.
    
    Requires a parameter file, use **harpy stitchparams** to generate one and adjust it for your study.
    The `--vcf-samples` option considers only the samples present in your input `--vcf` file rather than all
    the samples identified in `INPUT`. If providing additional STITCH arguments, they must be in quotes and 
    in R language style. Use single-quotes (string literals) if supplying an argument that requires quotes. For example:
    
    ```
    harpy ... -x 'switchModelIteration = 39, splitReadIterations = NA, reference_populations = c("CEU","GBR")'...
    ```
    """
    output_dir = output_dir.rstrip("/")
    workflowdir = f"{output_dir}/workflow"
    sdm = "conda" if conda else "conda apptainer"
    command = f'snakemake --rerun-incomplete --show-failed-logs --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/impute.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if snakemake is not None:
        command += snakemake

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    validate_input_by_ext(vcf, "--vcf", ["vcf", "bcf", "vcf.gz"])
    check_impute_params(parameters)
    bamlist, n = parse_alignment_inputs(inputs)
    validate_bam_RG(bamlist)
    samplenames = vcf_samplematch(vcf, bamlist, vcf_samples)
    biallelic, n_biallelic = biallelic_contigs(vcf, f"{workflowdir}")

    fetch_rule(workflowdir, "impute.smk")
    fetch_script(workflowdir, "stitch_impute.R")
    fetch_report(workflowdir, "impute.Rmd")
    fetch_report(workflowdir, "stitch_collate.Rmd")
    os.makedirs(f"{output_dir}/logs/snakemake", exist_ok = True)
    sm_log = snakemake_log(output_dir, "impute")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: impute\n")
        config.write(f"snakemake_log: {sm_log}\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samples_from_vcf: {vcf_samples}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skip_reports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  paramfile: {Path(parameters).resolve()}\n")
        config.write(f"  variantfile: {Path(vcf).resolve()}\n")
        config.write(f"  biallelic_contigs: {Path(biallelic).resolve()}\n")
        config.write("  alignments:\n")
        for i in bamlist:
            config.write(f"    - {i}\n")
    
    generate_conda_deps()
    if setup_only:
        sys.exit(0)

    start_text = Table(show_header=False,pad_edge=False, show_edge=False, padding = (0,0), box=box.SIMPLE)
    start_text.add_column("detail", justify="left", style="light_steel_blue", no_wrap=True)
    start_text.add_column("value", justify="left")
    start_text.add_row("Input VCF:", vcf)
    start_text.add_row("VCF Samples:", f"{len(samplenames)}")
    start_text.add_row("Alignment Files:", f"{n}")
    start_text.add_row("Parameter File:", parameters)
    start_text.add_row("Usable Contigs:", f"{n_biallelic}")
    start_text.add_row("Output Folder:", output_dir + "/")
    start_text.add_row("Workflow Log:", sm_log.replace(f"{output_dir}/", ""))
    launch_snakemake(command, "impute", start_text, output_dir, sm_log, quiet)
