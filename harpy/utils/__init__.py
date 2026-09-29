import click
from harpy import __version__
from harpy.common.lazy_group import LazyGroup, LazySpec

# Workflow rules call these many times, so each one is imported only when it is run
# (most pull in pysam, numpy, or polars). No help text is given, so listing them
# (e.g. `harpy-utils --help`) does import all of them.
COMMANDS = {
    "bx-stats-fq":       LazySpec("harpy.utils.bx_stats_fq:bx_stats_fq"),
    "bx-stats-sam":      LazySpec("harpy.utils.bx_stats_sam:bx_stats_sam"),
    "bx-to-end":         LazySpec("harpy.utils.bx_to_end:bx_to_end"),
    "check-bam":         LazySpec("harpy.utils.check_bam:check_bam"),
    "check-fastq":       LazySpec("harpy.utils.check_fastq:check_fastq"),
    "haplotag-acbd":     LazySpec("harpy.utils.haplotag_acbd:haplotag_acbd"),
    "infer-sv":          LazySpec("harpy.utils.infer_sv:infer_sv"),
    "known-adapters":    LazySpec("harpy.utils.known_adapters:known_adapters"),
    "molecule-coverage": LazySpec("harpy.utils.molecule_coverage:molecule_coverage"),
    "optical-dist-fq":   LazySpec("harpy.utils.optical_dist:optical_dist_fq"),
    "optical-dist-sam":  LazySpec("harpy.utils.optical_dist:optical_dist_sam"),
    "parse-phaseblocks": LazySpec("harpy.utils.parse_phaseblocks:parse_phaseblocks"),
    "plot-depth":        LazySpec("harpy.utils.plot_depth:plot_depth"),
    "process-notebook":  LazySpec("harpy.utils.process_notebook:process_notebook"),
    "rename-bam":        LazySpec("harpy.utils.rename_bam:rename_bam"),
}

@click.group(cls = LazyGroup, lazy_commands = COMMANDS, options_metavar='')
@click.version_option(__version__, prog_name="utils.hpy", hidden = True)
@click.help_option('--help', hidden = True)
def cli():
    "Utility scripts associated with Harpy"
