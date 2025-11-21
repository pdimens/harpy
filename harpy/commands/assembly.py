"""Perform a linked-read aware metassembly"""

from ssl import ALERT_DESCRIPTION_CLOSE_NOTIFY
import rich_click as click
import os
from harpy.common.cli_filetypes import HPCProfile, FASTQfile
from harpy.common.cli_types_generic import KParam, SnakemakeParams
from harpy.common.cli_types_params import SpadesParams, ArcsParams
from harpy.common.printing import workflow_info
from harpy.common.system_ops import container_ok
from harpy.common.workflow import Workflow
from harpy.validation.fastq import FASTQ

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/assembly")
# SPADES
@click.option('-k', '--kmer-length', panel = "Assembly Parameters", type = KParam(), show_default = True, default = "auto", help = 'K values to use for assembly (`odd` and `<128`)')
@click.option('-r', '--max-memory', panel = "Assembly Parameters",  type = click.IntRange(min = 1000), show_default = True, default = 10000, help = 'Maximum memory for spades to use, in megabytes')
@click.option('-x', '--extra-params', panel = "Assembly Parameters", type = SpadesParams(), help = 'Additional spades parameters, in quotes')
# TIGMINT/ARCS/LINKS
@click.option('-y', '--arcs-extra', panel = "Scaffolding Parameters", type = ArcsParams(), help = 'Additional ARCS parameters, in quotes (`option=arg` format)')
@click.option("-c", "--contig-length", panel = "Scaffolding Parameters", type = click.IntRange(min = 10), default = 500, show_default = True, help = "Minimum contig length")
@click.option("-n", "--links", panel = "Scaffolding Parameters", type = click.IntRange(min = 1), default = 5, show_default = True, help = "Minimum number of links to compute scaffold")
@click.option("-a", "--min-aligned", panel = "Scaffolding Parameters", type = click.IntRange(min = 1), default = 5, show_default = True, help = "Minimum aligned read pairs per barcode")
@click.option("-q", "--min-quality", panel = "Scaffolding Parameters", type = click.IntRange(0,40, clamp = True), default = 0, show_default = True, help = "Minimum mapping quality")
@click.option("-m", "--mismatch", panel = "Scaffolding Parameters", type = click.IntRange(min = 1), default = 5, show_default = True, help = "Maximum number of mismatches")
@click.option("-d", "--molecule-distance", panel = "Scaffolding Parameters", type = click.IntRange(min = 500), default = 50000, show_default = True, help = "Distance cutoff to split molecules (bp)")
@click.option("-l", "--molecule-length", panel = "Scaffolding Parameters", type = click.IntRange(min = 100), default = 2000, show_default = True, help = "Minimum molecule length (bp)")
@click.option("-i", "--seq-identity", panel = "Scaffolding Parameters", type = click.IntRange(0,100, clamp = True), default = 98, show_default = True, help = "Minimum sequence identity") 
@click.option("-s", "--span", panel = "Scaffolding Parameters", type = click.IntRange(min = 1), default = 20, show_default = True, help = "Minimum number of spanning molecules to be considered assembled")
# Other Options
@click.option('-o', '--output-dir', panel = "Workflow Options", type = click.Path(exists = False, resolve_path = True), default = "Assembly", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', panel = "Workflow Options", default = 4, show_default = True, type = click.IntRange(1, 999, clamp = True), help = 'Number of threads to use')
@click.option('-u', '--organism-type', panel = "Assembly Parameters", type = click.Choice(['prokaryote', 'eukaryote', 'fungus'], case_sensitive=False), default = "eukaryote", show_default=True, help = "Organism type for assembly report [`eukaryote`,`prokaryote`,`fungus`]")
@click.option('--clean', hidden = True, panel = "Workflow Options", type = str, help = 'Delete the log (`l`), .snakemake (`s`), and/or workflow (`w`) folders when done')
@click.option('--container', panel = "Workflow Options",  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--hpc', panel = "Workflow Options",  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', panel = "Workflow Options", default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` all output, `1` progress bar, `2` no output')
@click.option('--setup', panel = "Workflow Options",  is_flag = True, hidden = True, default = False, help = 'Setup the workflow and exit')
@click.option('--skip-reports', panel = "Workflow Options",  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', panel = "Workflow Options", type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('fastq_r1', required=True, type=FASTQfile(single=True), nargs=1)
@click.argument('fastq_r2', required=True, type=FASTQfile(single=True), nargs=1)
def assembly(fastq_r1, fastq_r2, kmer_length, max_memory, output_dir, extra_params,arcs_extra,contig_length,links,min_quality,min_aligned,mismatch,molecule_distance,molecule_length,seq_identity,span, organism_type, clean, container, threads, snakemake, quiet, hpc, setup, skip_reports):
    """
    Assemble linked reads into a genome

    The linked-read barcodes must be in `BX:Z` or `BC:Z` FASTQ header tags. If provided, values for `-k` must be
    separated by commas and without spaces (e.g. `-k 15,23,51`). It is strongly recommended to first deconvolve
    the input FASTQ files with `harpy deconvolve`.
    """
    workflow = Workflow("assembly", "assembly.smk", output_dir, container, clean, quiet)
    workflow.setup_snakemake(threads, hpc, snakemake)
    workflow.conda = ["assembly","qc"]

    ## checks and validations ##
    fastq = FASTQ([fastq_r1,fastq_r2], quiet= quiet > 0)

    workflow.reports["skip"] = skip_reports
    workflow.reports["organism-type"] = organism_type
    workflow.input(fastq.files[0], "fastq-r1")
    workflow.input(fastq.files[1], "fastq-r2")
    workflow.param('auto' if kmer_length == "auto" else ",".join(map(str,kmer_length)), "spades:k")
    workflow.param(max_memory, "spades:max-memory")
    if extra_params:
        workflow.param(extra_params, "spades:extra")
    workflow.param(min_quality, "tigmint:minimum-mapping-quality")
    workflow.param(mismatch, "tigmint:mismatch")
    workflow.param(molecule_distance, "tigmint:molecule-distance")
    workflow.param(molecule_length, "tigmint:molecule-length")
    workflow.param(span, "tigmint:span")
    workflow.param(min_aligned, "arcs:minimum-aligned-reads")
    workflow.param(contig_length, "arcs:minimum-contig-length")
    workflow.param(seq_identity, "arcs:minimum-sequence-identity")
    if arcs_extra:
        workflow.param(arcs_extra, "arcs:extra")
    workflow.param(links, "links:minimum-links")

    workflow.start_text = workflow_info(
        ("Kmer Length: ", workflow.parameters["spades"]["k"]),
        ("Output Folder:", os.path.relpath(output_dir) + "/")
    )

    workflow.initialize(setup)
