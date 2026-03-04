import re
import os
import sys
import click
import pysam

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.argument('platform', required = True, type=click.Choice(['haplotagging','stlfr','tellseq'], case_sensitive=False))
@click.argument('input', required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def check_fastq(platform, input):
    """
    File format validation for FASTQ file

    Parses a FASTQ file to check if any sequences don't conform to the SAM spec,
    whether BX:Z: is the last tag in the record, and the counts of: total reads,
    reads without BX:Z: tag, reads with incorrect barcode depending on the platform.
    """
    N_READS: int = 0
    NO_BX: int = 0
    BAD_BX: int = 0
    BAD_SAM_SPEC: int = 0
    BX_NOT_LAST: int = 0

    samspec = re.compile(r'[A-Z][A-Z]:[AifZHB]:')
    def check_samspec(fq_comment):
        splithead = fq_comment.split()
        for i in splithead:
            # if comments dont start with TAG:TYPE:, invalid SAM spec
            if not samspec.match(i):
                BAD_SAM_SPEC += 1
            # if the BX:Z: isn't at the end, add to BX_NOT_LAST
            if not splithead[-1].startswith('BX:Z'):
                BX_NOT_LAST += 1

    platform = platform.lower()
    if platform == "haplotagging":
        barcode = re.compile(r'A[0-9][0-9]C[0-9][0-9]B[0-9][0-9]D[0-9][0-9]')
        def check_read(fq_record):
            if 'BX:Z:' not in fq_record.comment:
                NO_BX += 1
                return
            if not barcode.search(fq_record.comment):
                BAD_BX += 1
            check_samspec(fq_record.comment)

    if platform == "stlfr":
        barcode = re.compile(r'#\d+_\d+_\d+$')
        def check_read(fq_record):
            if not barcode.search(fq_record.name):
                BAD_BX += 1
            check_samspec(fq_record.comment)

    if platform == "tellseq":
        barcode = re.compile(r'\:[ATCGN]+$')
        def check_read(fq_record):
            if not barcode.search(fq_record.name):
                BAD_BX += 1
            check_samspec(fq_record.comment)

    with pysam.FastxFile(input, persist=False) as fh:
        for entry in fh:
            N_READS += 1
            check_read(entry)

    values = [str(i) for i in [os.path.basename(input), N_READS, NO_BX, BAD_BX, BAD_SAM_SPEC, BX_NOT_LAST]]
    sys.stdout.write("\t".join(values) + "\n")
