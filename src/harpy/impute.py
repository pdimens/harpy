from .helperfunctions import fetch_file, generate_conda_deps, getnames, print_onstart
from .helperfunctions import vcfcheck, vcf_samplematch, biallelic_contigs
from .helperfunctions import check_impute_params, validate_bamfiles, parse_alignment_inputs
import rich_click as click
import subprocess
import sys
import os

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/impute/")
@click.option('-v', '--vcf', required = True, type=click.Path(exists=True, dir_okay=False),metavar = "File Path", help = 'Path to BCF/VCF file')
@click.option('-p', '--parameters', required = True, type=click.Path(exists=True, dir_okay=False), metavar = "File Path", help = 'STITCH parameter file (tab-delimited)')
@click.option('-x', '--extra-params', type = str, metavar = "String", help = 'Additional STITCH parameters, in quotes')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Use samples present in vcf file for imputation rather than those found the directory')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(min = 4, max_open = True), metavar = "Integer", help = 'Number of threads to use')
@click.option('-s', '--snakemake', type = str, metavar = "String", help = 'Additional Snakemake parameters, in quotes')
@click.option('-r', '--skipreports',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t generate any HTML reports')
@click.option('-q', '--quiet',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Don\'t show output text while running')
@click.option('--print-only',  is_flag = True, show_default = True, default = False, metavar = "Toggle", help = 'Print the generated snakemake command and exit')
@click.argument('input', required=True, type=click.Path(exists=True), nargs=-1)
def impute(input, parameters, threads, vcf, vcf_samples, extra_params, snakemake, skipreports, quiet, print_only):
    """
    Impute genotypes using variants and sequences
    
    Requires a parameter file, use **harpy stitchparams** to generate one and adjust it for your study.
    Use the `--vcf-samples` toggle to phase only the samples present in your input `--vcf` file rather than all
    the samples present in the `--directory`. Additional STITCH arguments, if provided, must be in quotes and 
    in R language style. The extra parameters will remain constant across different models.
    Use single-quotes (string literals) if supplying an argument that requires quotes. For example:
    
    ```harpy ... -x 'switchModelIteration = 39, splitReadIterations = NA, reference_populations = c("CEU","GBR")'```
    """
    command = f'snakemake --rerun-incomplete --nolock --use-conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append('Impute/workflow/impute.smk')
    command.append("--configfile")
    command.append("Impute/workflow/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)

    fetch_file("impute.smk", "Impute/workflow/")
    fetch_file("stitch_impute.R", "Impute/workflow/")
    for i in ["Impute", "StitchCollate"]:
        fetch_file(f"{i}.Rmd", "Impute/workflow/report/")
    ## validate inputs ##
    vcfcheck(vcf)
    ### check that samples in VCF match input directory
    sn = parse_alignment_inputs(input, "Impute/workflow/input")
    samplenames = vcf_samplematch(vcf, "Impute/workflow/input", vcf_samples)
    check_impute_params(parameters)
    validate_bamfiles("Impute/workflow/input", samplenames)

    # generate and store list of viable contigs (minimum of 2 biallelic SNPs)
    # doing it here so it doesn't have to run each time inside the workflow
    contigs = biallelic_contigs(vcf)

    with open("Impute/workflow/config.yml", "w") as config:
        config.write(f"seq_directory: Impute/workflow/input\n")
        config.write(f"samplenames: {samplenames}\n")
        config.write(f"variantfile: {vcf}\n")
        config.write(f"paramfile: {parameters}\n")
        config.write(f"contigs: {contigs}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    print_onstart(
        f"Input VCF: {vcf}\nSamples in VCF: {len(samplenames)}\nContigs Considered: {len(contigs)}\nOutput Directory: Impute/",
        "impute"
    )
    generate_conda_deps()
    _module = subprocess.run(command)
    sys.exit(_module.returncode)