#! /usr/bin/env python3

import sys
import rich_click as click

@click.group(no_args_is_help = True, deprecated=True)
def convert():
    """
    Convert between linked-read formats

    This module of Harpy has been deprecated and its function has been moved to the Djinn package,
    which should be provided with the standard conda-based Harpy installation.
    """
    sys.exit(1)