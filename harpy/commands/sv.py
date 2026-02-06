"""Harpy workflows to detect structural variants"""

import os
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTAfile, SAMfile
from harpy.common.cli_types_generic import ContigList, MultiInt, SnakemakeParams
from harpy.common.cli_types_params import LeviathanParams, NaibrParams
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow
from harpy.validation.fasta import FASTA
from harpy.validation.populations import Populations
from harpy.validation.xam import XAM

@click.group()
@click.help_option('--help', hidden = True)
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
@click.option('-x', '--extra-params', panel = "Parameters", type = LeviathanParams(), help = 'Additional leviathan parameters, in quotes')
@click.option('-i', '--iterations', panel = "Parameters", show_default = True, default=50, type = click.IntRange(min = 10), help = 'Number of iterations to perform through index (reduces memory)')
@click.option('-d', '--duplicates', panel = "Parameters", show_default = True, default=10, type = click.IntRange(min = 1), help = 'Consider SV of the same type as duplicates if their breakpoints are within this distance')
@click.option('-m', '--min-size', panel = "Parameters", type = click.IntRange(min = 10), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-s', '--sharing-thresholds', panel = "Parameters", type = MultiInt(3, minimum = 5, maximum = 100), default = "95,95,95", show_default=True, help = 'Percentile thresholds in the distributions of the number of shared barcodes for (small,medium,large) variants (no spaces)')
@click.option('-b', '--min-barcodes', panel = "Parameters", show_default = True, default=2, type = click.IntRange(min = 1), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "SV/leviathan", show_default=True,  help = 'Output directory name')
@click.option('-p', '--populations', panel = "Parameters", type=click.Path(exists = True, dir_okay=False, readable=True, resolve_path=True), help = 'File of `sample`_\\<TAB\\>_`population`')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--contigs', panel = "Workflow Options",  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.help_option('--help', hidden = True)
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def leviathan(inputs, output_dir, reference, min_size, min_barcodes, iterations, duplicates, sharing_thresholds, threads, populations, extra_params, snakemake, skip_reports, quiet, hpc, clean, container, contigs, setup):
    """
    Call structural variants using LEVIATHAN
    
    Provide the reference fasta followed by the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    Optionally specify `--populations` for population-pooled variant calling
    (**harpy template** can create that file). If you suspect Leviathan is missing certain variants
    you expect to find, try lowering `--sharing-thresholds`, _e.g._ `90,90,90`. The thresholds don't
    have to be the same across the different size classes.
    """
    workflow = Workflow("sv_leviathan", "sv_leviathan.smk", output_dir, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake)
    workflow.notebook_files = ["sv.ipynb"]
    workflow.conda = ["align", "variants"]

    ## checks and validations ##
    alignments = XAM(inputs)
    fasta = FASTA(reference)
    if contigs:
        fasta.match_contigs(contigs)     

    workflow.input(fasta.file, "reference")
    if populations:
        popfile = Populations(populations, alignments.files)
        popfile.copy_to_workflow(output_dir)
        workflow.input(popfile.file, "groupings")
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
        "Output Folder" : os.path.relpath(output_dir) + "/"
    }

    workflow.initialize(setup)

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/sv/naibr/")
@click.option('-x', '--extra-params', panel = "Parameters", type = NaibrParams(), help = 'Additional naibr parameters, in quotes')
@click.option('-b', '--min-barcodes', panel = "Parameters", show_default = True, default=2, type = click.IntRange(min = 1), help = 'Minimum number of barcode overlaps supporting candidate SV')
@click.option('-q', '--min-quality', panel = "Parameters", show_default = True, default=30, type = click.IntRange(min = 0, max = 40), help = 'Minimum mapping quality of reads to use')
@click.option('-m', '--min-size', panel = "Parameters", type = click.IntRange(min = 10), default = 1000, show_default=True, help = 'Minimum size of SV to detect')
@click.option('-d', '--molecule-distance', panel = "Parameters", default = 100000, show_default = True, type = click.IntRange(min = 100), help = 'Base-pair distance delineating separate molecules')
@click.option('-o', '--output-dir', panel = "Parameters", type = click.Path(exists = False, resolve_path = True), default = "SV/naibr", show_default=True,  help = 'Output directory name')
@click.option('-p', '--populations', panel = "Parameters", type=click.Path(exists = True, dir_okay=False, readable=True, resolve_path=True), help = 'File of `sample`_\\<TAB\\>_`population`')
@click.option('-t', '--threads',panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(4,999, clamp = True), help = 'Number of threads to use')
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--contigs', panel = "Workflow Options",  type = ContigList(), help = 'File or list of contigs to plot')
@click.option('--setup',  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.help_option('--help', hidden = True)
@click.argument('reference', type=FASTAfile(), required = True, nargs = 1)
@click.argument('inputs', required=True, type=SAMfile(), nargs=-1)
def naibr(inputs, output_dir, reference, min_size, min_barcodes, min_quality, threads, populations, molecule_distance, extra_params, snakemake, skip_reports, quiet, hpc, clean, container, contigs, setup):
    """
    Call structural variants using NAIBR
    
    Provide the reference fasta followed by the input alignment (`.bam`) files and/or directories at the end of the command as 
    individual files/folders, using shell wildcards (e.g. `data/drosophila*.bam`), or both.

    NAIBR requires **phased** bam files as input. This appears as the `HP` or `PS` tags
    in alignment records. If your bam files do not have either of these phasing tags
    (e.g. BWA/strobealign do not phase alignments), use `harpy phase bam` to do so.

    Optionally specify `--populations` for population-pooled variant calling (**harpy template** can create that file).
    """
    workflow = Workflow("sv_naibr", "sv_naibr.smk", output_dir, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake)
    workflow.notebook_files = ["sv.ipynb"]
    workflow.conda = ["variants"]

    ## checks and validations ##
    alignments = XAM(inputs, check_phase = True, quiet = quiet > 0)
    fasta =  FASTA(reference, quiet = quiet > 0)
    if contigs:
        fasta.match_contigs(contigs)

    workflow.notebooks["skip"] = skip_reports
    workflow.notebooks["plot-contigs"] = contigs if contigs else "default"
    workflow.input(fasta.file, "reference")
    if populations:
        popfile = Populations(populations, alignments.files)
        popfile.copy_to_workflow(output_dir)
        workflow.input(popfile.file, "groupings")
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
        "Output Folder" : os.path.relpath(output_dir) + "/"
    }    

    workflow.initialize(setup)

sv.add_command(leviathan)
sv.add_command(naibr)
