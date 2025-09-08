"""Harpy haplotype phasing workflow"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, SAMfile, VCFfile, FASTAfile
from harpy.validation.fasta import FASTA
from harpy.validation.sam import SAM
from harpy.validation.vcf import VCF
from harpy.common.cli_types_generic import ContigList, SnakemakeParams
from harpy.common.cli_types_params import HapCutParams
from harpy.common.printing import workflow_info
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow

docstring = {
        "harpy phase": [
        {
            "name": "Parameters",
            "options": ["--extra-params", "--reference", "--min-map-quality", "--min-base-quality", "--molecule-distance", "--prune-threshold", "--unlinked"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--contigs", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads",  "--vcf-samples", "--help"],
            "panel_styles": {"border_style": "dim"}
        }
    ]
}

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/phase")
@click.option('-x', '--extra-params', type = HapCutParams(), help = 'Additional HapCut2 parameters, in quotes')
@click.option('-r', '--reference', type=FASTAfile(), help = 'Path to reference genome if wanting to also extract reads spanning indels')
@click.option('-q', '--min-map-quality', default = 20, show_default = True, type = click.IntRange(0, 40, clamp = True), help = 'Minimum mapping quality for phasing')
@click.option('-m', '--min-base-quality', default = 13, show_default = True, type = click.IntRange(0, 100, clamp = True), help = 'Minimum base quality for phasing')
@click.option('-d', '--molecule-distance', default = 100000, show_default = True, type = click.IntRange(min = 100), help = 'Distance cutoff to split molecules (bp)')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Phase", show_default=True,  help = 'Output directory name')
@click.option('-p', '--prune-threshold', default = 30, show_default = True, type = click.IntRange(0,100, clamp = True), help = 'PHRED-scale threshold (%) for pruning low-confidence SNPs (larger prunes more.)')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(2, 999, clamp = True), help = 'Number of threads to use')
@click.option('-U','--unlinked', is_flag = True, default = False, help = "Treat input data as not linked reads")
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--contigs',  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` unified progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--vcf-samples',  is_flag = True, show_default = True, default = False, help = 'Use samples present in vcf file for phasing rather than those found in the inputs')
@click.argument('vcf', required = True, type = VCFfile(), nargs = 1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def phase(vcf, inputs, output_dir, threads, unlinked, min_map_quality, min_base_quality, molecule_distance, prune_threshold, vcf_samples, reference, snakemake, extra_params, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Phase SNPs into haplotypes

    Provide the vcf file followed by the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/myotis*.bam`), or both.
    
    Presence and type of linked-read data is auto-detected, but you may choose to omit barcode
    information with `-U`. Use `--vcf-samples` to phase only the samples present in your input
    `VCF` file rather than all the samples present in the `INPUT` alignments.
    """
    workflow = Workflow("phase", "phase.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["hapcut.qmd"]
    workflow.conda = ["phase", "report"]

    ## checks and validations ##
    alignments = SAM(inputs, detect_bc= not unlinked)
    vcffile = VCF(vcf, workflow.workflow_directory)
    vcffile.match_samples(alignments.files, vcf_samples)

    if reference:
        fasta = FASTA(reference)
    if contigs:
        vcffile.match_contigs(contigs)

    workflow.config = {
        "workflow" : workflow.name,
        "phasing" : {
            "prune" : prune_threshold,
            "min_map_quality": min_map_quality,
            "min_base_quality": min_base_quality
        },
        "linkedreads": {
            "type" : alignments.lr_type,
            "distance_threshold" : molecule_distance,
        },
        **({'extra': extra_params} if extra_params else {}),
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "reports" : {
            "skip": skip_reports,
            **({'plot_contigs': contigs} if contigs else {'plot_contigs': "default"}),
        },
        "inputs" : {
            "vcf" : {
                "file": vcffile.file,
                "prioritize_samples" : vcf_samples
            },
            **({'reference': fasta.file} if reference else {}),
            "alignments" : alignments.files
        }
    }

    workflow.start_text = workflow_info(
        ("Input VCF:", os.path.basename(vcffile.file)),
        ("Samples:", min(len(vcffile.samples), alignments.count)),
        ("Barcode Type:", alignments.lr_type),
        ("Phase Indels:", "yes" if reference else "no"),
        ("Reference:", os.path.basename(reference)) if reference else None,
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)
