"""Harpy workflows to detect structural variants"""

from harpy.validation.vcf import VCF
from harpy.common.cli_filetypes import VCFfile
from harpy.common.cli_filetypes import FASTQfile
import os

import rich_click as click

from harpy.common.cli_filetypes import FASTAfile, FASTQfile, HPCProfile, SAMfile
from harpy.common.cli_params import (
    ContigList,
    LeviathanParams,
    MultiInt,
    NaibrParams,
    SnakemakeParams,
)
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow
from harpy.validation.fasta import FASTA
from harpy.validation.fastq import FASTQ
from harpy.validation.populations import Populations
from harpy.validation.xam import XAM

@click.group()
@click.help_option('--help', hidden = True)
def sv():
    """
    Call inversions, deletions, and duplications from alignments

    |   caller  | inversions | duplications | deletions | breakends |
    |:----------|:----------:|:------------:|:---------:|:---------:|
    | leviathan |      ✔     |     ✔        |     ✔     |      ✔    |
    | naibr     |      ✔     |     ✔        |     ✔     |     🗙     |

    Provide the subcommand `leviathan` or `naibr` to get more information on using
    those variant callers. Once structural variants are identified, use `genotype` to
    genotype them using `svjedi-tag`.
    """

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog= "Documentation: https://pdimens.github.io/harpy/workflows/sv/leviathan/")
@click.option('-x', '--extra-params', panel = "Parameters", type = LeviathanParams(), help = 'Additional leviathan parameters, in quotes')
@click.option('-i', '--iterations', panel = "Parameters", show_default = True, default=50, type = click.IntRange(min = 10), help = 'Number of iterations to perform through index (reduces memory)')
@click.option('-d', '--duplicates', panel = "Parameters", show_default = True, default=10, type = click.IntRange(min = 1), help = 'Consider SV of the same type as duplicates if their breakpoints are within this distance')
@click.option('-m', '--min-size', panel = "Parameters", type = click.IntRange(min = 10), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-s', '--sharing-thresholds', panel = "Parameters", type = MultiInt(3, minimum = 5, maximum = 100), default = "95,95,95", show_default=True, help = 'Percentile thresholds in the distributions of the number of shared barcodes for (small,medium,large) variants (no spaces)')
@click.option('-b', '--min-barcodes', panel = "Parameters", show_default = True, default=2, type = click.IntRange(min = 1), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-p', '--populations', panel = "Parameters", type=click.Path(exists = True, dir_okay=False, readable=True, resolve_path=True), help = 'File of `sample`_\\<TAB\\>_`population`')
@click.option('-O', '--output', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "SV/leviathan", show_default=True,  help = 'Output directory name')
@click.option('-@', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-T', '--no-temp', hidden = True, panel = "Workflow Options", is_flag = True, default = False, help = 'Don\'t delete temporary files')
@click.option('-C', '--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('-N', '--setup',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('-H', '--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('-Q', '--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('-R', '--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('-S', '--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--contigs', panel = "Workflow Options",  type = ContigList(), help = 'File or list of contigs to plot')
@click.help_option('--help', hidden = True)
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def leviathan(inputs, output, reference, min_size, min_barcodes, iterations, duplicates, sharing_thresholds, threads, populations, extra_params, snakemake, skip_reports, quiet, hpc, clean, container, contigs, setup, no_temp):
    """
    Call structural variants using LEVIATHAN

    Provide the reference fasta followed by the input alignment (`.bam`) files and/or directories at the end of the command as
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    Optionally specify `--populations` for population-pooled variant calling
    (**harpy template** can create that file). If you suspect Leviathan is missing certain variants
    you expect to find, try lowering `--sharing-thresholds`, _e.g._ `90,90,90`. The thresholds don't
    have to be the same across the different size classes.
    """
    workflow = Workflow("sv_leviathan", "sv_leviathan.smk", output, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake, no_temp)
    workflow.notebook_files = ["sv.ipynb"]
    workflow.conda = ["align", "variants"]

    ## checks and validations ##
    alignments = XAM(inputs, detect_bc = True, nonlinked_ok = False, quiet = quiet)
    fasta = FASTA(reference, quiet)
    if contigs:
        fasta.match_contigs(contigs)

    workflow.input(fasta.file, "reference")
    if populations:
        popfile = Populations(populations, alignments.files, quiet)
        popfile.copy_to_workflow(output)
        workflow.input(popfile.file, "groupings:source")
        workflow.input("workflow/sample.groups", "groupings:processed")
    workflow.input(alignments.files, "alignments")

    workflow.notebooks["skip"] = skip_reports
    workflow.notebooks["plot-contigs"] = contigs if contigs else "default"
    workflow.param(min_barcodes, "min-barcodes")
    workflow.param(min_size, "min-size")
    workflow.param(iterations, "iterations")
    workflow.param(sharing_thresholds[0], "variant-thresholds:small")
    workflow.param(sharing_thresholds[1], "variant-thresholds:medium")
    workflow.param(sharing_thresholds[2], "variant-thresholds:large")
    workflow.param(duplicates, "variant-thresholds:duplicates")
    if extra_params:
        workflow.param(extra_params, "extra")

    workflow.info = {
        "Samples" : alignments.count,
        "Reference" : os.path.basename(reference),
        "Sample Pooling" : os.path.basename(populations) if populations else "no",
        "Output Folder" : os.path.relpath(output) + "/"
    }

    workflow.initialize(setup)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/sv/naibr/")
@click.option('-x', '--extra-params', panel = "Parameters", type = NaibrParams(), help = 'Additional naibr parameters, in quotes')
@click.option('-b', '--min-barcodes', panel = "Parameters", show_default = True, default=2, type = click.IntRange(min = 1), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-q', '--min-quality', panel = "Parameters", show_default = True, default=30, type = click.IntRange(min = 0, max = 40), help = 'Minimum mapping quality of reads to use')
@click.option('-m', '--min-size', panel = "Parameters", type = click.IntRange(min = 10), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-d', '--molecule-distance', panel = "Parameters", default = 100000, show_default = True, type = click.IntRange(min = 100), help = 'Base-pair distance delineating separate molecules')
@click.option('-p', '--populations', panel = "Parameters", type=click.Path(exists = True, dir_okay=False, readable=True, resolve_path=True), help = 'File of `sample`_\\<TAB\\>_`population`')
@click.option('-O', '--output', panel = "Parameters", type = click.Path(exists = False, resolve_path = True), default = "SV/naibr", show_default=True,  help = 'Output directory name')
@click.option('-@', '--threads',panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-T', '--no-temp', hidden = True, panel = "Workflow Options", is_flag = True, default = False, help = 'Don\'t delete temporary files')
@click.option('-C', '--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('-N', '--setup',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('-H', '--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('-Q', '--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('-R', '--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('-S', '--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--contigs', panel = "Workflow Options",  type = ContigList(), help = 'File or list of contigs to plot')
@click.help_option('--help', hidden = True)
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def naibr(inputs, output, reference, min_size, min_barcodes, min_quality, threads, populations, molecule_distance, extra_params, snakemake, skip_reports, quiet, hpc, clean, container, contigs, setup, no_temp):
    """
    Call structural variants using NAIBR

    Provide the reference fasta followed by the input alignment (`.bam`) files and/or directories at the end of the command as
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    NAIBR requires **phased** bam files as input. This appears as the `HP` or `PS` tags
    in alignment records. If your bam files do not have either of these phasing tags
    (e.g. BWA/strobealign do not phase alignments), use `harpy phase bam` to do so.

    Optionally specify `--populations` for population-pooled variant calling (**harpy template** can create that file).
    """
    workflow = Workflow("sv_naibr", "sv_naibr.smk", output, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake, no_temp)
    workflow.notebook_files = ["sv.ipynb"]
    workflow.conda = ["variants"]

    ## checks and validations ##
    alignments = XAM(inputs, detect_bc = True, nonlinked_ok = False, check_phase = True, quiet = quiet)
    fasta =  FASTA(reference, quiet = quiet)
    if contigs:
        fasta.match_contigs(contigs)

    workflow.notebooks["skip"] = skip_reports
    workflow.notebooks["plot-contigs"] = contigs if contigs else "default"
    workflow.input(fasta.file, "reference")
    if populations:
        popfile = Populations(populations, alignments.files, quiet)
        popfile.copy_to_workflow(output)
        workflow.input(popfile.file, "groupings:source")
        workflow.input("workflow/sample.groups", "groupings:processed")
    workflow.input(alignments.files, "alignments")
    workflow.param(min_barcodes, "min-barcodes")
    workflow.param(min_quality, "min-map-quality")
    workflow.param(min_size, "min-size")
    workflow.param(molecule_distance, "molecule-distance")
    if extra_params:
        workflow.param(extra_params, "extra")

    workflow.info = {
        "Samples" : alignments.count,
        "Reference" : os.path.basename(reference),
        "Sample Pooling" : os.path.basename(populations) if populations else "no",
        "Output Folder" : os.path.relpath(output) + "/"
    }

    workflow.initialize(setup)


@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/sv/svjeditag/")
@click.option("-r", "--region-size", panel = "Parameters", default=10000, show_default = True, type = click.IntRange(min = 50), help = "size of barcode scanning window (bp)")
@click.option("-i", "--inaccuracy", panel = "Parameters", default=0, show_default = True, type = click.IntRange(min = 0), help = "breakpoint inaccuracy tolerance (bp)")
@click.option("-m", "--min-diff", panel = "Parameters",  default=20, show_default = True, type = click.IntRange(min = 1), help = "minimum difference between the two best genotype likelihoods")
@click.option("-e", "--prob-error", panel = "Parameters", default=[0.2, 0.1, 0.02, 0.008], show_default=True, type = click.FloatRange(min=0.0, max = 1.0), nargs = 4, help = "genotype likelihood probability error according to inversion size [small medium large xl]")
#@click.option('-g', '--gfa', panel = "Parameters", show_default = True, type = click.Path(exists = True, dir_okay=False, resolve_path=True), help = 'GFA graph file previously generated by SVJedi-Tag')
@click.option('-O', '--output', panel = "Parameters", type = click.Path(exists = False, resolve_path = True), default = "SV/genotype", show_default=True,  help = 'Output directory name')
@click.option('-@', '--threads',panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('-T', '--no-temp', hidden = True, panel = "Workflow Options", is_flag = True, default = False, help = 'Don\'t delete temporary files')
@click.option('-C', '--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('-N', '--setup',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('-H', '--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('-Q', '--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('-S', '--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.help_option('--help', hidden = True)
@click.argument('vcf', required = True, type = VCFfile(), nargs = 1)
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=FASTQfile(), nargs=-1)
def genotype(vcf, reference, inputs, output, region_size, inaccuracy, min_diff, prob_error, threads, snakemake, quiet, hpc, clean, container, setup, no_temp):
    """
    Genotype called variants using SVJedi-Tag

    SVJedi-Tag does not identify structural variants-- it genotypes previously called ones, such as the output of LEVIATHAN or NAIBR.
    Provide the vcf of SVs and reference fasta, followed by one or more paired-end FASTQ files at the end of the command.
    If setting custom probability error rates, you must put in 4 values (small med large xl) e.g., `--e 0.2 0.1 0.02 0.008`, where the values correspond to different SV sizes,
    given by: small: <25kb, medium 25kb-50kb, large 50kb-100kb, xl >=100kb. 
    """
    workflow = Workflow("sv_genotype", "sv_genotype.smk", output, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake, no_temp)
    #workflow.notebook_files = ["sv.ipynb"]
    workflow.conda = ["variants"]

    ## checks and validations ##
    fastq = FASTQ(inputs, detect_bc = False, nonlinked_ok = False, quiet = quiet)
    fastq.has_bx_tag()
    fasta =  FASTA(reference, quiet = quiet)
    _vcf = VCF(vcf, workdir=workflow.workflow_directory, quiet = quiet, threads = threads)
    _vcf.is_sv_fmt()
    _vcf.get_contigs()
    fasta.match_contigs_vcf(_vcf.contigs, _vcf.file)

    #workflow.notebooks["skip"] = skip_reports
    workflow.input(fasta.file, "reference")
    workflow.input(_vcf.file, "vcf")
    workflow.input(fastq.files, "fastq")
    #if gfa:
    #    workflow.input(gfa, "gfa")
    workflow.param(region_size, "region-size")
    workflow.param(inaccuracy, "inaccuracy")
    workflow.param(min_diff, "likelihood-min-diff")
    workflow.param(list(prob_error), "likelihood-prob-error")

    workflow.info = {
        "Reference" : os.path.basename(reference),
        "Variants" : os.path.basename(vcf),
        "Search Window" : region_size,
        "Output Folder" : os.path.relpath(output) + "/"
    }

    workflow.initialize(setup)


sv.add_command(leviathan)
sv.add_command(naibr)
sv.add_command(genotype)
