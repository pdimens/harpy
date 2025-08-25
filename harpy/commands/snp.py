"""Harpy workflows to call SNP variants"""

import os
import shutil
from pathlib import Path
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTAfile, PopulationFile, SAMfile
from harpy.common.cli_types_generic import SnakemakeParams, SNPRegion
from harpy.common.cli_types_params import MpileupParams, FreebayesParams
from harpy.common.system_ops import container_ok
from harpy.common.printing import workflow_info
from harpy.common.workflow import Workflow
from harpy.validation.fasta import FASTA
from harpy.validation.populations import Populations
from harpy.validation.sam import SAM

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


@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/snp")
@click.option('-x', '--extra-params', type = FreebayesParams(), help = 'Additional freebayes parameters, in quotes')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path= True), default = "SNP/freebayes", show_default=True,  help = 'Output directory name')
@click.option('-n', '--ploidy', default = 2, show_default = True, type=click.IntRange(min=1), help = 'Ploidy of samples')
@click.option('-p', '--populations', type=PopulationFile(), help = 'File of `sample`_\\<TAB\\>_`population`')
@click.option('-r', '--regions', type=SNPRegion(), default=50000000, show_default=True, help = "Regions where to call variants")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` unified progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
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
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["bcftools_stats.qmd"]
    workflow.conda = ["r", "variants"]

    ## checks and validations ##
    alignments = SAM(inputs)
    fasta = FASTA(reference)
    fasta.validate_region(regions)

    region = Path(os.path.join(workflow.workflow_directory, "regions.bed")).resolve().as_posix()
    if isinstance(regions, int):
        os.system(f"make_windows -m 1 -w {regions} {reference} > {region}")
    elif os.path.exists(regions):
        shutil.copy2(regions, region)
    else:
        region = regions
    if populations:
        popfile = Populations(populations, alignments.files)

    workflow.config = {
        "workflow" : workflow.name,
        "ploidy" : ploidy,
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "reports" : {"skip": skip_reports},
        "inputs" : {
            "reference" : fasta.file,
            "regions" : region,
            **({'groupings': popfile.file} if populations else {}),
            "alignments" : alignments.files
        }
    }

    workflow.start_text = workflow_info(
        ("Samples:", alignments.count),
        ("Sample Groups:", os.path.basename(populations)) if populations else None,
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/snp")
@click.option('-x', '--extra-params', type = MpileupParams(), help = 'Additional mpileup parameters, in quotes')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path=True), default = "SNP/mpileup", show_default=True,  help = 'Output directory name')
@click.option('-n', '--ploidy', default = 2, show_default = True, type=click.IntRange(1, 2), help = 'Ploidy of samples')
@click.option('-p', '--populations', type=PopulationFile(), help = 'File of `sample`\\<TAB\\>`population`')
@click.option('-r', '--regions', type=SNPRegion(), default=50000000, show_default=True, help = "Regions where to call variants")
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--quiet', show_default = True, default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` unified progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def mpileup(reference, inputs, output_dir, regions, threads, populations, ploidy, extra_params, snakemake, skip_reports, quiet, hpc, container, setup_only):
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
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["bcftools_stats.qmd"]
    workflow.conda = ["r"]

    ## checks and validations ##
    alignments = SAM(inputs)
    fasta = FASTA(reference)
    fasta.validate_region(regions)

    region = Path(os.path.join(workflow.workflow_directory, "regions.bed")).resolve().as_posix()
    if isinstance(regions, int):
        os.system(f"make_windows -m 1 -w {regions} {reference} > {region}")
    elif os.path.exists(regions):
        shutil.copy2(regions, region)
    else:
        region = regions
    if populations:
        popfile = Populations(populations, alignments.files)

    workflow.config = {
        "workflow" : workflow.name,
        "ploidy" : ploidy,
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "reports" : {
            "skip": skip_reports
        },
        "inputs" : {
            "reference" : fasta.file,
            "regions" : region,
            **({'groupings': popfile.file} if populations else {}),
            "alignments" : alignments.files
        }
    }

    workflow.start_text = workflow_info(
        ("Samples:", alignments.count),
        ("Sample Groups:", os.path.basename(populations)) if populations else None,
        ("Reference:", os.path.basename(reference)),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

snp.add_command(mpileup)
snp.add_command(freebayes)
