"""Harpy imputation workflow"""

import os
import sys
import subprocess
import rich_click as click
from pathlib import Path
from .conda_deps import generate_conda_deps
from .helperfunctions import fetch_rule, fetch_report, fetch_script
from .fileparsers import parse_alignment_inputs, biallelic_contigs
from .printfunctions import print_onstart
from .validations import validate_input_by_ext, vcf_samplematch, check_impute_params, validate_bam_RG

docstring = {
        "harpy impute": [
        {
            "name": "Parameters",
            "options": ["--vcf", "--parameters", "--extra-params", "--vcf-samples"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--hpc", "--conda", "--snakemake", "--quiet", "--help"],
        },
    ]
}

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/impute/")
@click.option('-v', '--vcf', required = True, type=click.Path(exists=True, dir_okay=False),metavar = "File Path", help = 'Path to BCF/VCF file')
@click.option('-p', '--parameters', required = True, type=click.Path(exists=True, dir_okay=False), help = 'STITCH parameter file (tab-delimited)')
@click.option('-x', '--extra-params', type = str, help = 'Additional STITCH parameters, in quotes')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for imputation rather than those found the inputs')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), help = 'Number of threads to use')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, help = 'Don\'t show output text while running')
@click.option('-o', '--output-dir', type = str, default = "Impute", show_default=True, help = 'Name of output directory')
@click.option('--snakemake', type = str, help = 'Additional Snakemake parameters, in quotes')
@click.option('--hpc',  type = click.Path(exists = True, file_okay = False), help = 'Config dir for automatic HPC submission')
@click.option('--conda',  is_flag = True, default = False, help = 'Use conda/mamba instead of container')
@click.option('--skipreports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--print-only',  is_flag = True, hidden = True, default = False, help = 'Print the generated snakemake command and exit')
@click.argument('inputs', required=True, type=click.Path(exists=True), nargs=-1)
def impute(inputs, output_dir, parameters, threads, vcf, vcf_samples, extra_params, snakemake, skipreports, quiet, hpc, conda, print_only):
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
    command = f'snakemake --rerun-incomplete --rerun-triggers input mtime params --nolock --software-deployment-method {sdm} --conda-prefix ./.snakemake/conda --cores {threads} --directory . '
    command += f"--snakefile {workflowdir}/impute.smk "
    command += f"--configfile {workflowdir}/config.yaml "
    if hpc:
        command += f"--workflow-profile {hpc} "
    if quiet:
        command += "--quiet all "
    if snakemake is not None:
        command += snakemake
    if print_only:
        click.echo(command)
        sys.exit(0)

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    validate_input_by_ext(vcf, "--vcf", ["vcf", "bcf", "vcf.gz"])
    check_impute_params(parameters)
    bamlist, n = parse_alignment_inputs(inputs)
    validate_bam_RG(bamlist)
    samplenames = vcf_samplematch(vcf, bamlist, vcf_samples)
    biallelic, n_biallelic = biallelic_contigs(vcf, f"{workflowdir}")

    fetch_rule(workflowdir, "impute.smk")
    fetch_script(workflowdir, "stitch_impute.R")
    for i in ["Impute", "StitchCollate"]:
        fetch_report(workflowdir, f"{i}.Rmd")

    with open(f"{workflowdir}/config.yaml", "w", encoding="utf-8") as config:
        config.write("workflow: impute\n")
        config.write(f"output_directory: {output_dir}\n")
        config.write(f"samples_from_vcf: {vcf_samples}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {command}\n")
        config.write("inputs:\n")
        config.write(f"  paramfile: {Path(parameters).resolve()}\n")
        config.write(f"  variantfile: {Path(vcf).resolve()}\n")
        config.write(f"  biallelic_contigs: {Path(biallelic).resolve()}\n")
        config.write("  alignments:\n")
        for i in bamlist:
            config.write(f"    - {i}\n")
    
    print_onstart(
        f"Input VCF: {vcf}\nSamples in VCF: {len(samplenames)}\nAlignments Provided: {n}\nContigs with ≥2 Biallelic SNPs: {n_biallelic}\nOutput Directory: {output_dir}/",
        "impute"
    )
    generate_conda_deps()
    _module = subprocess.run(command.split())
    sys.exit(_module.returncode)
