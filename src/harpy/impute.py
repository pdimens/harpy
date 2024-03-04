from .helperfunctions import fetch_file, generate_conda_deps, biallelic_contigs
from .fileparsers import parse_alignment_inputs
from .printfunctions import print_onstart
from .validations import vcfcheck, vcf_samplematch, check_impute_params, validate_bamfiles
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
    
    Provide the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.
    
    Requires a parameter file, use **harpy stitchparams** to generate one and adjust it for your study.
    Use the `--vcf-samples` toggle to phase only the samples present in your input `--vcf` file rather than all
    the samples present in the `--directory`. Additional STITCH arguments, if provided, must be in quotes and 
    in R language style. The extra parameters will remain constant across different models.
    Use single-quotes (string literals) if supplying an argument that requires quotes. For example:
    
    ```
    harpy ... -x 'switchModelIteration = 39, splitReadIterations = NA, reference_populations = c("CEU","GBR")'...
    ```
    """
    workflowdir = "Impute/workflow"
    command = f'snakemake --rerun-incomplete --nolock --software-deployment-method conda --conda-prefix ./.snakemake/conda --cores {threads} --directory .'.split()
    command.append('--snakefile')
    command.append(f'{workflowdir}/impute.smk')
    command.append("--configfile")
    command.append(f"{workflowdir}/config.yml")
    if quiet:
        command.append("--quiet")
        command.append("all")
    if snakemake is not None:
        [command.append(i) for i in snakemake.split()]
    call_SM = " ".join(command)
    if print_only:
        click.echo(call_SM)
        exit(0)

    os.makedirs(f"{workflowdir}/", exist_ok= True)
    sn = parse_alignment_inputs(input, f"{workflowdir}/input")
    samplenames = vcf_samplematch(vcf, f"{workflowdir}/input", vcf_samples)
    validate_bamfiles(f"{workflowdir}/input", samplenames)
    ## validate inputs ##
    vcfcheck(vcf)
    check_impute_params(parameters)
    fetch_file("impute.smk", f"{workflowdir}/")
    fetch_file("stitch_impute.R", f"{workflowdir}/")
    for i in ["Impute", "StitchCollate"]:
        fetch_file(f"{i}.Rmd", f"{workflowdir}/report/")

    # generate and store list of viable contigs (minimum of 2 biallelic SNPs)
    # doing it here so it doesn't have to run each time inside the workflow
    contigs = biallelic_contigs(vcf)
    with open(f"{workflowdir}/config.yml", "w") as config:
        config.write(f"seq_directory: {workflowdir}/input\n")
        config.write(f"samplenames: {samplenames}\n")
        config.write(f"variantfile: {vcf}\n")
        config.write(f"paramfile: {parameters}\n")
        config.write(f"contigs: {contigs}\n")
        if extra_params is not None:
            config.write(f"extra: {extra_params}\n")
        config.write(f"skipreports: {skipreports}\n")
        config.write(f"workflow_call: {call_SM}\n")

    generate_conda_deps()
    print_onstart(
        f"Input VCF: {vcf}\nSamples in VCF: {len(samplenames)}\nAlignments Provided: {len(sn)}\nContigs Considered: {len(contigs)}\nOutput Directory: Impute/",
        "impute"
    )
    _module = subprocess.run(command)
    sys.exit(_module.returncode)