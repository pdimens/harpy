"""Downsample a fastq/bam file by barcodes"""

import sys
import rich_click as click

@click.command(no_args_is_help = True, deprecated=True)
def downsample():
    """
    Downsample data by barcode

    This module of Harpy has been deprecated and its function has been moved to the Djinn package,
    which should be provided with the standard conda-based Harpy installation.
    """
    sys.exit(1)