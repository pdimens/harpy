"""Harpy workflows to call SNP variants"""

import os
import sys
from pathlib import Path
import rich_click as click
from .common.cli_types_generic import HPCProfile, InputFile, SnakemakeParams, SNPRegion
from .common.cli_types_params import MpileupParams, FreebayesParams
from .common.parsers import parse_alignment_inputs
from .common.printing import workflow_info
from .common.validations import check_fasta, validate_bam_RG, validate_popfile, validate_popsamples, validate_regions
from .common.workflow import Workflow

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def snp():
    """
    Call SNPs and small indels from alignments
    
    Provide an additional subcommand `mpileup` or `freebayes` to get more information on using
    those variant callers. They are both robust variant callers, but `freebayes` is recommended when ploidy
    is greater than **2**.
    """

module_docstring = {
    "harpy snp": [
        {
            "name": "Commands",
            "commands": ["freebayes", "mpileup"],
            "panel_styles": {"border_style": "blue"}
        }
    ]
}

docstring = {
    "harpy snp mpileup": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--ploidy", "--populations", "--regions"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ],
    "harpy snp freebayes": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--ploidy", "--populations", "--regions"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
} |  module_docstring


@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/snp")
@click.option('-x', '--extra-params', type = FreebayesParams(), help = 'Additional freebayes parameters, in quotes')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path= True), default = "SNP/freebayes", show_default=True,  help = 'Output directory name')
@click.option('-n', '--ploidy', default = 2, show_default = True, type=click.IntRange(min=1), help = 'Ploidy of samples')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False, readable=True, resolve_path=True), help = 'File of `sample`_\\<TAB\\>_`population`')
@click.option('-r', '--regions', type=SNPRegion(), default=50000000, show_default=True, help = "Regions where to call variants")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=InputFile("fasta", gzip_ok = True), required = True, nargs = 1)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def freebayes(reference, inputs, output_dir, threads, populations, ploidy, regions, extra_params, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Call variants using freebayes
    
    Provide the reference fasta followed by the input alignment (`.bam`) files and/or directories
    at the end of the command as individual files/folders, using shell wildcards
    (e.g. `data/jellyfish*.bam`), or both.
    
    The `--regions` option specifies what genomic regions to call variants
    with. If a 1-indexed BED file is provided, variant calling will be parallelized
    over those regions. If a single region is provided in the format `chrom:start-end`, only
    that region will be called. If an integer is provided (default), then Harpy will
    call variants in parallel for intervals of that size across the entire reference genome.

    Optionally specify `--populations` for population-aware variant calling (**harpy template** can create that file).
    """
    workflow = Workflow("snp_freebayes", "snp_freebayes.smk", output_dir, quiet)
    workflow.reports = ["bcftools_stats.qmd"]
    workflow.conda = ["r", "variants"]

    ## checks and validations ##
    bamlist, n = parse_alignment_inputs(inputs, "INPUTS")
    validate_bam_RG(bamlist, threads, quiet)
    check_fasta(reference)
    validate_regions(regions, reference)
    region = Path(f"{workflowdir}/regions.bed").resolve().as_posix()
    if isinstance(regions, int):
        os.system(f"make_windows.py -m 1 -w {regions} {reference} > {region}")
    elif os.path.exists(regions):
        shutil.copy2(regions, region)
    else:
        region = regions
    if populations:
        # check for delimeter and formatting
        validate_popfile(populations)
        # check that samplenames and populations line up
        validate_popsamples(bamlist, populations,quiet)

    ## workflow setup ##
    workflow.setup_snakemake(
        "conda" if not container else "conda apptainer",
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    workflow.config = {
        "workflow" : workflow.name,
        "ploidy" : ploidy,
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "absolute": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "reports" : {"skip": skip_reports},
        "inputs" : {
            "reference" : reference,
            "regions" : region,
            **({'groupings': populations} if populations else {}),
            "alignments" : bamlist
        }
    }

    workflow.start_text = workflow_info(
        ("Samples:", n),
        ("Sample Groups:", os.path.basename(populations)) if populations else None,
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize()
    if not setup_only:
        workflow.launch("workflow/snp.freebayes.summary")

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/snp")
@click.option('-x', '--extra-params', type = MpileupParams(), help = 'Additional mpileup parameters, in quotes')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path=True), default = "SNP/mpileup", show_default=True,  help = 'Output directory name')
@click.option('-n', '--ploidy', default = 2, show_default = True, type=click.IntRange(1, 2), help = 'Ploidy of samples')
@click.option('-p', '--populations', type=click.Path(exists = True, dir_okay=False, readable=True, resolve_path=True), help = 'File of `sample`\\<TAB\\>`population`')
@click.option('-r', '--regions', type=SNPRegion(), default=50000000, show_default=True, help = "Regions where to call variants")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=InputFile("fasta", gzip_ok = True), required = True, nargs = 1)
@click.argument('inputs', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=-1)
def mpileup(inputs, output_dir, regions, reference, threads, populations, ploidy, extra_params, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Call variants from using bcftools mpileup
    
    Provide the reference fasta followed by the input alignment (`.bam`) files and/or directories
    at the end of the command as individual files/folders, using shell wildcards
    (e.g. `data/scarab*.bam`), or both.
    
    The `--regions` option specifies what genomic regions to call variants
    with. If a 1-indexed BED file is provided, variant calling will be parallelized
    over those regions. If a single region is provided in the format `chrom:start-end`, only
    that region will be called. If an integer is provided (default), then Harpy will
    call variants in parallel for intervals of that size across the entire reference genome.

    Optionally specify `--populations` for population-aware variant calling (**harpy template** can create that file).
    """
    workflow = Workflow("snp_mpileup", "snp_mpileup.smk", output_dir, quiet)
    workflow.reports = ["bcftools_stats.qmd"]
    workflow.conda = ["r"]

    ## checks and validations ##
    bamlist, n = parse_alignment_inputs(inputs, "INPUTS")
    validate_bam_RG(bamlist, threads, quiet)
    check_fasta(reference)
    validate_regions(regions, reference)
    region = Path(f"{workflowdir}/regions.bed").resolve().as_posix()
    if isinstance(regions, int):
        os.system(f"make_windows.py -m 1 -w {regions} {reference} > {region}")
    elif os.path.isfile(regions):
        shutil.copy2(regions, region)
    else:
        region = regions
    if populations:
        validate_popfile(populations)
        # check that samplenames and populations line up
        validate_popsamples(bamlist, populations, quiet)

    ## workflow setup ##
    workflow.setup_snakemake(
        "conda" if not container else "conda apptainer",
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    workflow.config = {
        "workflow" : workflow.name,
        "ploidy" : ploidy,
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "absolute": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "reports" : {
            "skip": skip_reports
        },
        "inputs" : {
            "reference" : reference,
            "regions" : region,
            **({'groupings': populations} if populations else {}),
            "alignments" : bamlist
        }
    }

    workflow.start_text = workflow_info(
        ("Samples:", n),
        ("Sample Groups:", os.path.basename(populations)) if populations else None,
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize()
    if not setup_only:
        workflow.launch("workflow/snp.mpileup.summary")

snp.add_command(mpileup)
snp.add_command(freebayes)
