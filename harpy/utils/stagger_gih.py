import click
import pysam
import sys
import io
from itertools import zip_longest
import signal
if hasattr(signal, "SIGPIPE"):
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

pad  = ["TTTTTTT", "CCCCCC", "GGGGG", "AAAA", "TTT", "CC", "GG", ""]
qpad = ["I" * len(i) for i in pad]

@click.help_option('--help', hidden=True)
@click.command(hidden = True, no_args_is_help=True, epilog="Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.argument('fastq', nargs=2, required=True, type=click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.argument("info", type=click.File("r"), default="-", required=False)
def stagger_gih(fastq, info):
    """
    Given an info file from cutadapt, adds any necessary staggers to the barcodes
    in the input FASTQ file. Writes to stdout. INTENDED FOR INTERNAL USE in harpy preprocess gih.
    """
    t = 0
    discarded = 0

    # buffer stdout to reduce syscalls
    stdout = io.open(sys.stdout.fileno(), mode="wb", buffering=2**19, closefd=False)
    write = stdout.write

    with (
        pysam.FastxFile(fastq[0], persist=False) as fq_F,
        pysam.FastxFile(fastq[1], persist=False) as fq_R,
    ):
        for cutadapt, fq1, fq2 in zip_longest(info, fq_F, fq_R):
            t += 1
            if cutadapt and fq1:
                fields = cutadapt.rstrip('\n').split('\t')
                CA_name = fields[0].partition(' ')[0]
                col2 = int(fields[1])
            else:
                break

            if col2 == -1:
                discarded += 1
                continue

            if CA_name != fq1.name:
                sys.stderr.write(f"Error: Read name mismatch at record {t}. Expected:\n{CA_name}\nbut found:\n{fq1.name}\n")
                sys.exit(1)

            col3 = int(fields[2])
            plen = 7 if (col3 < 51 or col3 > 58) else (col3 - 51)

            p = pad[plen]
            q = qpad[plen]

            if plen == 6:
                seq  = p + fq1.sequence[1:]
                qual = q + fq1.quality[1:]
            else:
                seq  = p + fq1.sequence
                qual = q + fq1.quality

            # write FASTQ records directly as bytes, bypassing object construction
            write(f"@{fq1.name} {fq1.comment}\n{seq}\n+\n{qual}\n{fq2}\n".encode())

    stdout.flush()
    sys.stderr.write(f"Reads discarded because of missing ME sequence: {discarded}\n")