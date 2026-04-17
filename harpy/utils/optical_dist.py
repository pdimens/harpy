import click
from pysam import AlignmentFile, FastxFile
import sys

seqCodes = {
    "HW": 100,
    "M" : 100, 
    "K" : 100,
    "N" : 100,
    "E" : 100,
    "ST-E" : 100,
    "A" : 2500,
    "V" : 2500,
    "H" : 2500,
    "LH" : 2500
}

@click.command(hidden = True, no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.argument("bam", required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def optical_dist_sam(bam):
    '''
    Read the first record of a BAM file and print the optical duplication distance parameter (100 or 2500)
    based on the instrument code of the sequence name. INTERNAL USE ONLY.
    '''
    with AlignmentFile(bam, require_index=False) as inBam:
        for record in inBam.fetch(until_eof=True):
            prefix = record.query_name.partition(":")[0]
            for i in seqCodes:
                if prefix.startswith(i):
                    sys.stdout.write(f"{seqCodes[i]}\n")
                    sys.exit(0)
            sys.stdout.write(f"{seqCodes[i]}\n")
            sys.exit(0)

@click.command(hidden = True, no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.argument("fastq", nargs = -1, required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def optical_dist_fq(fastq):
    '''
    Read the first record of a FASTQ file and print the optical duplication distance parameter (100 or 2500)
    based on the instrument code of the sequence name. INTERNAL USE ONLY.
    '''
    with FastxFile(fastq[0], persist=False) as fq:
        for record in fq:
            _id = record.name
            for i in seqCodes:
                if _id.startswith(i):
                    sys.stdout.write(f"{seqCodes[i]}\n")
                    sys.exit(0)
            sys.stdout.write(f"{seqCodes[i]}\n")
            sys.exit(0)