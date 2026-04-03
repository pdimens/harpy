import click
from harpy import __version__

from .bx_stats_sam import bx_stats_sam
from .bx_stats_fq import bx_stats_fq
from .bx_to_end import bx_to_end
from .check_bam import check_bam
from .check_fastq import check_fastq
from .haplotag_acbd import haplotag_acbd
from .infer_sv import infer_sv
from .molecule_coverage import molecule_coverage
from .parse_phaseblocks import parse_phaseblocks
from .plot_depth import plot_depth
from .process_notebook import process_notebook
from .rename_bam import rename_bam
from .optical_dist import optical_dist_sam, optical_dist_fq

@click.group(options_metavar='')
@click.version_option(__version__, prog_name="utils.hpy", hidden = True)
@click.help_option('--help', hidden = True)
def cli():
    "Utility scripts associated with Harpy"

cli.add_command(bx_stats_sam)
cli.add_command(bx_to_end)
cli.add_command(check_bam)
cli.add_command(check_fastq)
cli.add_command(bx_stats_fq)
cli.add_command(haplotag_acbd)
cli.add_command(infer_sv)
cli.add_command(molecule_coverage)
cli.add_command(optical_dist_sam)
cli.add_command(optical_dist_fq)
cli.add_command(parse_phaseblocks)
cli.add_command(plot_depth)
cli.add_command(process_notebook)
cli.add_command(rename_bam)