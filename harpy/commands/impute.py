"""Harpy imputation workflow"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, SAMfile, VCFfile, ImputeStrategy
from harpy.common.cli_params import SnakemakeParams, StitchParams
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow
from harpy.validation.impute_parameters import ImputeParams
from harpy.validation.xam import XAM
from harpy.validation.vcf import VCF
##@click.option('-r', '--region', panel = "Parameters", type = str, help = 'Specific region to impute')

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/impute/")
@click.option('-x', '--extra-params', panel = "Parameters", type = StitchParams(), help = 'Additional STITCH parameters, in quotes')
@click.option('-b', '--buffer', panel = "Parameters", default = 100000, show_default = True, type = click.IntRange(min=1, clamp = True), help = 'Base pairs to consider on each side of genomic `region` or `window`')
@click.option('-g', '--grid-size', panel = "Parameters", hidden = True, show_default = True, default = 1, type = click.IntRange(min = 1), help = 'Perform imputation in windows of a specific size, instead of per-SNP (default)')
@click.option('-O', '--output', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Impute", show_default=True,  help = 'Output directory name')
@click.option('-s', '--strategy', panel = "Parameters", type = ImputeStrategy(), default = "window:1000000", help = 'Imputation strategy (see above)')
@click.option('-@', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(2,999, clamp = True), help = 'Number of threads to use')
@click.option('-T', '--no-temp', hidden = True, panel = "Workflow Options", is_flag = True, default = False, help = 'Don\'t delete temporary files')
@click.option('-C', '--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('-N', '--setup',  panel = "Workflow Options", is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('-H', '--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('-Q', '--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('-R', '--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('-S', '--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('-V', '--vcf-samples', panel = "Parameters",  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for imputation rather than those found in the inputs')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.help_option('--help', hidden = True)
@click.argument('parameters', required = True, type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True), nargs=1)
@click.argument('vcf', required = True, type = VCFfile(), nargs=1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def impute(parameters, vcf, inputs, output, strategy, buffer, grid_size, threads, vcf_samples, extra_params, snakemake, skip_reports, quiet, hpc, clean, container, setup, no_temp):
    """
    Impute variant genotypes from alignments
    
    Provide the parameter file followed by the input VCF and the input alignment files/directories (`.bam`) at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.
    
    Use `harpy template impute` to generate a parameter file. Enclose any additional STITCH arguments
    in quotes and with the `--option=value` format, without spaces (e.g. `"--switchModelIteration=39"`).
    
    # Imputation Strategies (--strategy)
    ## Buffered genomic windows (default)
    Imputation can be very memory-intensive, so the default strategy is imputation in 1Mb windows with a 100kb buffer
    (100kb before and after window). This takes the format `window:size`, e.g., `window:1000000`
    ## A specifc chromosomal region
    Use the format: `contig:start-end` to impute only a single genomic region, e.g., `ch1:1-500000`.
    ## All chromosomes in their entirety
    The default approach in previous harpy versions and is triggered by the word `all` and ignores the `--buffer` value.
    This can be very memory intensive and is no longer recommended.
    """
    workflow = Workflow("impute", "impute.smk", output, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake, no_temp)
    workflow.notebook_files = ["impute.ipynb", "stitch_collate.ipynb"]
    workflow.conda = ["impute"]

    ## checks and validations ##
    params = ImputeParams(parameters, quiet)
    alignments = XAM(inputs, quiet = quiet)
    vcffile = VCF(vcf, workflow.workflow_directory, quiet)
    vcffile.find_biallelic_contigs()
    vcffile.match_samples(alignments.files, vcf_samples)
    region = strategy.lower() != "all" and not strategy.lower().startswith("window:")
    window = None
    if region:
        cntg,start,end = vcffile.validate_region(strategy)
        region = f"{cntg}:{start}-{end}"
    elif strategy.lower().startswith("window:"):
        window = int(strategy.split(":")[1])

    workflow.notebooks["skip"] = skip_reports
    workflow.input(params.file, "parameters")
    workflow.input(vcffile.file, "vcf")
    if not region:
        workflow.input(vcffile.biallelic_file, "biallelic-contigs")
    workflow.input(alignments.files, "alignments")
    workflow.param(grid_size, "grid-size")
    if region:
        workflow.param(region, "region")
    elif window:
        workflow.param(window, "window-size")
    if strategy != "all":
        workflow.param(buffer, "buffer")
    if extra_params:
        workflow.param(extra_params, "extra")
    workflow.param(params.parameters, "stitch")

    workflow.info = {
        "Input VCF" : os.path.basename(vcf),
        "Samples": min(len(vcffile.samples), alignments.count),
        "Parameter File" : os.path.basename(parameters),
        **({'Contigs': f"{len(vcffile.contigs)} [dim](with ≥ 5 biallelic SNPs)"} if not region else {"Target Region" : region}),
        **({'Window Size': window} if window else {}),
        **({'Buffer Size': buffer} if strategy != "all" else {}),
        'Output Folder' : os.path.relpath(output) + "/"
    }

    workflow.initialize(setup)
