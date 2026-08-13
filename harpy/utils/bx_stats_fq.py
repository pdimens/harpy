import sys

import click
from pysam import FastxFile

from harpy.common.parsers import Haplotagging, Stlfr, Tellseq

@click.command(no_args_is_help = True)
@click.argument('platform', required = True, type=click.Choice(['haplotagging','stlfr','tellseq'], case_sensitive=False))
@click.argument('input', required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def bx_stats_fq(platform, input):
    """
    Parses a FASTQ file to count: total sequences, total number of linked-read barcodes,
    number of valid barcodes, number of invalid BX tags, and a count of positional
    barcode invalidations (e.g. A00, _0_, N)
    """
    if platform == "haplotagging":
        parser = Haplotagging()
    elif platform == "stlfr":
        parser = Stlfr()
    else:
        parser = Tellseq()

    N_READS = 0
    N_BX = 0
    N_VALID = 0

    with FastxFile(input, persist = False) as fh:
        for entry in fh:
            N_READS += 1
            query = parser.getBX_fq(entry)
            if query:
                N_BX += 1
                barcode = query.group(1)
                is_invalid = parser.process_invalid(barcode)
                if not is_invalid:
                    N_VALID += 1

    sys.stdout.write(f"Reads\t{N_READS}\n")
    sys.stdout.write(f"Barcodes\t{N_BX}\n")
    sys.stdout.write(f"Barcodes_Valid\t{N_VALID}\n")
    sys.stdout.write(f"Barcodes_Invalid\t{N_BX - N_VALID}\n")
    for bx,count in parser.invalidDict.items():
        sys.stdout.write(f"{bx}\t{count}\n")
