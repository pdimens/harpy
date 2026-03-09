import click
import pysam
import sys
from itertools import zip_longest
import signal
if hasattr(signal, "SIGPIPE"):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def format_rec(rec, padlen) -> None:
    '''Process a record and give it padding for the stagger, writing to stdout'''
    if padlen == 6:
        seq = pad[padlen] + rec.sequence[1:]
        qual = qpad[padlen] + rec.quality[1:]
    else:
        seq = pad[padlen] + rec.sequence
        qual = qpad[padlen] + rec.quality
    
    sys.stdout.write(
        str(pysam.FastxRecord(rec.name, rec.comment, seq, qual)) + "\n"
    )

pad = ["TTTTTTT", "CCCCCC", "GGGGG", "AAAA", "TTT", "CC", "GG", ""]
qpad = ["I" * len(i) for i in pad]

@click.help_option('--help', hidden = True)
@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.argument('fastq', nargs=2, required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.argument("info", type=click.File("r"), default="-", required=False)
def stagger_gih(fastq, info):
    """
    Given an info file from cutadapt and a processed info_summary file, adds any necessary staggers
    to the barcodes in the input FASTQ file. Writes to stdout. INTENDED FOR INTERNAL USE in harpy preprocess gih.
    """
    t = 0
    discarded = 0
    with (
        pysam.FastxFile(fastq[0], persist=False) as fq_F,
        pysam.FastxFile(fastq[1], persist=False) as fq_R
    ):
        for cutadapt, fq1, fq2 in zip_longest(info, fq_F, fq_R):
            t += 1
            if cutadapt and fq1:
                fields = cutadapt.rstrip('\n').split('\t')
                CA_name = fields[0].split(' ')[0]
                col2 = int(fields[1])
            else:
                # end when there are no more lines in R1
                break
            if col2 == -1:
                discarded += 1
                continue

            if CA_name != fq1.name:
                print(f"Error: Read name mismatch at record {t}. Expected:\n{CA_name}\nbut found:\n{fq1.name}")
                sys.exit(1)

            col3 = int(fields[2])
            plen = 7 if (col3 < 51 or col3 > 58) else (col3 - 51)

            format_rec(fq1, plen)
            sys.stdout.write(str(fq2) + '\n')

    sys.stderr.write(f"Reads discarded because of missing ME sequence: {discarded}\n")
