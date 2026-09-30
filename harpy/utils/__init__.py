import click
from harpy import __version__
from harpy.common.lazy_group import LazyGroup, LazySpec

# Workflow rules call these many times, so each one is imported only when it is run
# (most pull in pysam, numpy, or polars). The help text (first paragraph of each command's docstring)
# and hidden flag are repeated here so `--help` and shell completion can list them without importing
# anything. test/cli/check_lazy_cli.py checks that this stays in sync with the commands.
COMMANDS = {
    "bx-stats-fq":       LazySpec("harpy.utils.bx_stats_fq:bx_stats_fq", "Parses a FASTQ file to count: total sequences, total number of linked-read barcodes, number of valid barcodes, number of invalid BX tags, and a count of positional barcode invalidations (e.g. A00, _0_, N)"),
    "bx-stats-sam":      LazySpec("harpy.utils.bx_stats_sam:bx_stats_sam", "Linked-read metrics from alignment files"),
    "bx-to-end":         LazySpec("harpy.utils.bx_to_end:bx_to_end", "Move BX:Z tag to the end of records"),
    "check-bam":         LazySpec("harpy.utils.check_bam:check_bam", "File format validation for SAM/BAM file"),
    "check-fastq":       LazySpec("harpy.utils.check_fastq:check_fastq", "File format validation for FASTQ file"),
    "haplotag-acbd":     LazySpec("harpy.utils.haplotag_acbd:haplotag_acbd", "Generates the BC_{ABCD}.txt files necessary to demultiplex Gen I haplotagging barcodes"),
    "infer-sv":          LazySpec("harpy.utils.infer_sv:infer_sv", "Infer variant types from NAIBR bedpe output"),
    "known-adapters":    LazySpec("harpy.utils.known_adapters:known_adapters", "INTERNAL USE- Writes a fasta file of common illumina adapters that might appear in the GIH Nextera prep, the one used by Cornell GIH for tagmentation, and the ME sequence for use in QC adapter trimming. Writes to stdout. Most of these were derived from fastp (https://github.com/OpenGene/fastp)"),
    "molecule-coverage": LazySpec("harpy.utils.molecule_coverage:molecule_coverage", "Calculate molecule coverage from a barcode stats file"),
    "optical-dist-fq":   LazySpec("harpy.utils.optical_dist:optical_dist_fq", "Read the first record of a FASTQ file and print the optical duplication distance parameter (100 or 2500) based on the instrument code of the sequence name. INTERNAL USE ONLY.", hidden = True),
    "optical-dist-sam":  LazySpec("harpy.utils.optical_dist:optical_dist_sam", "Read the first record of a BAM file and print the optical duplication distance parameter (100 or 2500) based on the instrument code of the sequence name. INTERNAL USE ONLY.", hidden = True),
    "parse-phaseblocks": LazySpec("harpy.utils.parse_phaseblocks:parse_phaseblocks", "Summarize a HapCut2 phase block file"),
    "plot-depth":        LazySpec("harpy.utils.plot_depth:plot_depth", "Plot histograms of alignment and/or molecule depths"),
    "process-notebook":  LazySpec("harpy.utils.process_notebook:process_notebook", "Replace placeholder text in jupyter notebooks"),
    "rename-bam":        LazySpec("harpy.utils.rename_bam:rename_bam", "Rename a SAM/BAM file and modify the @RG tag"),
}

@click.group(cls = LazyGroup, lazy_commands = COMMANDS, options_metavar='')
@click.version_option(__version__, prog_name="utils.hpy", hidden = True)
@click.help_option('--help', hidden = True)
def cli():
    "Utility scripts associated with Harpy"
