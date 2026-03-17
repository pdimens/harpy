import re
import sys
import os
import click
import pysam

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.argument('platform', required = True, type=click.Choice(['10x','haplotagging','stlfr','tellseq'], case_sensitive=False))
@click.argument('input', required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def check_bam(platform, input):
    """
    File format validation for SAM/BAM file
    
    Specific to linked-read data. Checks if the sample name matches the RG tag,
    whether BX:Z: is the last tag in the record, and the counts of: total alignments,
    alignments with an MI:i: tag, alignments without BX:Z: tag, incorrect BX:Z: tag. Writes to stdout.
    """
    if platform == "haplotagging":
        bc_pattern = re.compile(r'^A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]')
    elif platform == "stlfr":
        bc_pattern = re.compile(r'^\d+_\d+_\d+')
    else:
        bc_pattern = re.compile(r'^[ATCGN]+')

    bam_in = input
    filename = os.path.basename(bam_in)
    bam_pattern = re.compile(r"\.[bB][aA][mM]$", flags = re.IGNORECASE)
    corename = re.sub(bam_pattern, "", filename)

    N_READS   = 0
    NO_BX      = 0
    BAD_BX     = 0
    BX_NOT_LAST = 0
    NO_MI      = 0
    NAME_MISMATCH = 0

    with pysam.AlignmentFile(bam_in, require_index=False) as alnfile:
        if alnfile.header.get("RG")[0]['ID'] != corename:
            NAME_MISMATCH += 1

        for record in alnfile.fetch(until_eof=True):
            N_READS += 1
            try:
                bx = record.get_tag("BX")
            except KeyError:
                NO_BX += 1
                continue
            # do a regex search to find proper barcode pattern in the BX
            if not re.search(bc_pattern, bx):
                # malformed BX tag
                BAD_BX += 1
            # do a search to see if BX:Z: tag is last tag in record
            if record.get_tags()[-1][0] != 'BX':
                BX_NOT_LAST += 1
            try:
                record.get_tag("MI")
            except KeyError:
                NO_MI += 1

    values = [str(i) for i in [filename, N_READS, NAME_MISMATCH, NO_MI, NO_BX, BX_NOT_LAST, BAD_BX]]
    sys.stdout.write("\t".join(values) + "\n")
