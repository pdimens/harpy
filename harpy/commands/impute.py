"""Harpy imputation workflow"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, SAMfile, VCFfile
from harpy.common.cli_types_generic import SnakemakeParams
from harpy.common.cli_types_params import StitchParams
from harpy.common.system_ops import container_ok
from harpy.common.printing import workflow_info
from harpy.common.workflow import Workflow
from harpy.validation.impute_parameters import ImputeParams
from harpy.validation.sam import SAM
from harpy.validation.vcf import VCF

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/impute/")
@click.option('-x', '--extra-params', panel = "Parameters", type = StitchParams(), help = 'Additional STITCH parameters, in quotes')
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Impute", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(2,999, clamp = True), help = 'Number of threads to use')
@click.option('-r', '--region', panel = "Parameters", type = str, help = 'Specific region to impute')
@click.option('-g', '--grid-size', panel = "Parameters", show_default = True, default = 1, type = click.IntRange(min = 1), help = 'Perform imputation in windows of a specific size, instead of per-SNP (default)')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only',  panel = "Workflow Options", is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--vcf-samples', panel = "Parameters",  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for imputation rather than those found in the inputs')
@click.argument('parameters', required = True, type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True), nargs=1)
@click.argument('vcf', required = True, type = VCFfile(), nargs=1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def impute(parameters, vcf, inputs, output_dir, region, grid_size, threads, vcf_samples, extra_params, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Impute variant genotypes from alignments
    
    Provide the parameter file followed by the input VCF and the input alignment files/directories (`.bam`) at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.
    
    Use `harpy template` to generate one and adjust it for your study. Set a `--grid-size` (in base pairs)
    to significantly reduce computation time and memory usage at the cost of minor accuracy loss.
    The `--vcf-samples` option considers only the samples present in your input `VCF` file rather than all
    the samples identified in `INPUTS`. Use `--region` to only impute a specific genomic region, given as
    `contig:start-end-buffer`, otherwise all contigs will be imputed. If providing additional STITCH arguments, they
    must be in quotes and in the `--option=value` format, without spaces (e.g. `"--switchModelIteration=39"`).
    """
    workflow = Workflow("impute", "impute.smk", output_dir, container, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake)
    workflow.reports = ["impute.qmd", "stitch_collate.qmd"]
    workflow.conda = ["report", "stitch"]

    ## checks and validations ##
    params = ImputeParams(parameters)
    alignments = SAM(inputs)
    vcffile = VCF(vcf, workflow.workflow_directory)
    vcffile.find_biallelic_contigs()
    vcffile.match_samples(alignments.files, vcf_samples)
    if region:
        vcffile.validate_region(region)

    workflow.inputs = {
        "parameters" : params.file,
        "vcf" : vcffile.file,
        **({"biallelic_contigs" : vcffile.biallelic_file} if not region else {}), 
        "alignments" : alignments.files
    }
    workflow.config = {
        "workflow" : workflow.name,
        **({'stitch_extra': extra_params} if extra_params else {}),
        **({'region': region} if region else {}),
        "reports" : {"skip": skip_reports},
        "grid_size": grid_size,
        "stitch_parameters" : params.parameters,
    }

    workflow.start_text = workflow_info(
        ("Input VCF:", os.path.basename(vcf)),
        ("Samples:", min(len(vcffile.samples), alignments.count)),
        ("Parameter File:", os.path.basename(parameters)),
        ("Contigs:", f"{len(vcffile.biallelic_contigs)} [dim](with at least 5 biallelic SNPs)") if not region else ("Target Region:", region),
        ("Output Folder:", os.path.relpath(output_dir) + "/")
    )

    workflow.initialize(setup_only)
