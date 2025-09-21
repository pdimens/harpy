"""Harpy workflows to detect structural variants"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTAfile, SAMfile, VCFfile
from harpy.common.cli_types_generic import ContigList, MultiInt, PANEL_OPTIONS, SnakemakeParams
from harpy.common.cli_types_params import LeviathanParams, NaibrParams
from harpy.common.printing import workflow_info
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow
from harpy.validation.fasta import FASTA
from harpy.validation.populations import Populations
from harpy.validation.sam import SAM
from harpy.validation.vcf import VCF

@click.group(context_settings={"help_option_names" : []})
def sv():
    """
    Call inversions, deletions, and duplications from alignments
 
    | caller | inversions | duplications | deletions | breakends |
    |:-------|:----------:|:------------:|:---------:|:---------:|
    | leviathan |      ✔  |     ✔        |     ✔     |      ✔    |
    | naibr     |      ✔  |     ✔        |     ✔     |     🗙    |

    Provide the subcommand `leviathan` or `naibr` to get more information on using
    those variant callers. NAIBR tends to call variants better, but requires more user preprocessing.
    """

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/sv/leviathan/")
@click.rich_config(PANEL_OPTIONS)
@click.option_panel("Parameters", panel_styles = {"border_style": "blue"})
@click.option_panel("Workflow Options", options = ["--help"],   panel_styles = {"border_style": "blue"})
@click.option('-x', '--extra-params', panel = "Parameters", type = LeviathanParams(), help = 'Additional leviathan parameters, in quotes')
@click.option('-i', '--iterations', panel = "Parameters", show_default = True, default=50, type = click.IntRange(min = 10), help = 'Number of iterations to perform through index (reduces memory)')
@click.option('-d', '--duplicates', panel = "Parameters", show_default = True, default=10, type = click.IntRange(min = 1), help = 'Consider SV of the same type as duplicates if their breakpoints are within this distance')
@click.option('-m', '--min-size', panel = "Parameters", type = click.IntRange(min = 10), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-s', '--sharing-thresholds', panel = "Parameters", type = MultiInt(3, minimum = 5, maximum = 100), default = "95,95,95", show_default=True, help = 'Percentile thresholds in the distributions of the number of shared barcodes for (small,medium,large) variants (no spaces)')
@click.option('-b', '--min-barcodes', panel = "Parameters", show_default = True, default=2, type = click.IntRange(min = 1), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "SV/leviathan", show_default=True,  help = 'Output directory name')
@click.option('-p', '--populations', panel = "Parameters", type=click.Path(exists = True, dir_okay=False, readable=True, resolve_path=True), help = 'File of `sample`_\\<TAB\\>_`population`')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--contigs', panel = "Workflow Options",  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def leviathan(inputs, output_dir, reference, min_size, min_barcodes, iterations, duplicates, sharing_thresholds, threads, populations, extra_params, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Call structural variants using LEVIATHAN
    
    Provide the reference fasta followed by the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    Optionally specify `--populations` for population-pooled variant calling
    (**harpy template** can create that file). If you suspect Leviathan is missing certain variants
    you expect to find, try lowering `--sharing-thresholds`, _e.g._ `90,90,90`. The thresholds don't
    have to be the same across the different size classes.
    """
    vcaller = "sv_leviathan" if not populations else "sv_leviathan_pop"
    workflow = Workflow("sv_leviathan", f"{vcaller}.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["leviathan.qmd"]
    if populations:
        workflow.reports.append("leviathan_pop.qmd")
    workflow.conda = ["align", "report", "variants"]

    ## checks and validations ##
    alignments = SAM(inputs)
    fasta = FASTA(reference)
    if contigs:
        fasta.match_contigs(contigs)
    if populations:
        popfile = Populations(populations, alignments.files)

    workflow.config = {
        "workflow" : workflow.name,
        "min_barcodes" : min_barcodes,
        "min_size" : min_size,
        "iterations" : iterations,
        "variant_thresholds": {
            "small" : sharing_thresholds[0],
            "medium" : sharing_thresholds[1],
            "large" : sharing_thresholds[2],
            "duplicates": duplicates
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
            "reference" : fasta.file,
            **({'groupings': popfile.file} if populations else {}),
            "alignments" : alignments.files
        }
    }

    workflow.start_text = workflow_info(
        ("Samples:", alignments.count),
        ("Reference:", os.path.basename(reference)),
        ("Sample Pooling:", os.path.basename(populations) if populations else "no"),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/sv/naibr/")
@click.rich_config(PANEL_OPTIONS)
@click.option_panel("Parameters", panel_styles = {"border_style": "blue"})
@click.option_panel("Workflow Options", options = ["--help"], panel_styles = {"border_style": "blue"})
@click.option('-x', '--extra-params', panel = "Parameters", type = NaibrParams(), help = 'Additional naibr parameters, in quotes')
@click.option('-b', '--min-barcodes', panel = "Parameters", show_default = True, default=2, type = click.IntRange(min = 1), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-q', '--min-quality', panel = "Parameters", show_default = True, default=30, type = click.IntRange(min = 0, max = 40), help = 'Minimum mapping quality of reads to use')
@click.option('-m', '--min-size', panel = "Parameters", type = click.IntRange(min = 10), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-d', '--molecule-distance', panel = "Parameters", default = 100000, show_default = True, type = click.IntRange(min = 100), help = 'Base-pair distance delineating separate molecules')
@click.option('-o', '--output-dir', panel = "Parameters", type = click.Path(exists = False, resolve_path = True), default = "SV/naibr", show_default=True,  help = 'Output directory name')
@click.option('-p', '--populations', panel = "Parameters", type=click.Path(exists = True, dir_okay=False, readable=True, resolve_path=True), help = 'File of `sample`_\\<TAB\\>_`population`')
@click.option('-t', '--threads',panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-v', '--vcf', panel = "Parameters", type=VCFfile(),  help = 'Path to phased bcf/vcf file')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--contigs', panel = "Workflow Options",  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup-only',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def naibr(inputs, output_dir, reference, vcf, min_size, min_barcodes, min_quality, threads, populations, molecule_distance, extra_params, snakemake, skip_reports, quiet, hpc, container, contigs, setup_only):
    """
    Call structural variants using NAIBR
    
    Provide the reference fasta followed by the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    NAIBR requires **phased** bam files as input. This appears as the `HP` or `PS` tags
    in alignment records. If your bam files do not have either of these phasing tags
    (e.g. BWA/strobealign do not phase alignments), then provide a **phased** `--vcf` file such
     as that created by `harpy phase` and Harpy will use [whatshap haplotag](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag)
    to phase your input bam files prior to calling variants with NAIBR.

    Optionally specify `--populations` for population-pooled variant calling (**harpy template** can create that file).
    """
    vcaller = "sv_naibr" if not populations else "sv_naibr_pop"
    vcaller += "_phase" if vcf else ""
    workflow = Workflow("sv_naibr", f"{vcaller}.smk", output_dir, quiet)
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.reports = ["naibr.qmd"]
    if populations:
        workflow.reports.append("naibr_pop.qmd")
    workflow.conda = ["phase", "report", "variants"]

    ## checks and validations ##
    alignments = SAM(inputs)
    fasta =  FASTA(reference)
    if contigs:
        fasta.match_contigs(contigs)
    if populations:
        popfile = Populations(populations, alignments.files)
    if vcf:
        vcffile = VCF(vcf, workflow.workflow_directory)
        vcffile.check_phase()

    workflow.config = {
        "workflow" : workflow.name,
        "min_barcodes" : min_barcodes,
        "min_quality" : min_quality,
        "min_size" : min_size,
        "molecule_distance" : molecule_distance,
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
            **({'reference': fasta.file} if reference else {}),
            **({'vcf': vcffile.file} if vcf else {}),
            **({'groupings': popfile.file} if populations else {}),
            "alignments" : alignments.files
        }
    }

    workflow.start_text = workflow_info(
        ("Samples:", alignments.count),
        ("Reference:", os.path.basename(reference)),
        ("Sample Pooling:", os.path.basename(populations) if populations else "no"),
        ("Perform Phasing:", "yes" if vcf else "no"),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

sv.add_command(leviathan)
sv.add_command(naibr)
