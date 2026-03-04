import click
import pysam
import sys
from itertools import zip_longest

def needs_stagger(filename):
    '''
    Using a summary generated from the cutadapt-generated INFO file, determine
    if the barcodes need a stagger inserted to make their lengths consistent for
    demuxing.
    '''
    counts = {}
    total_count = 0
    in_col3 = False

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line == 'col3=startpost':
                in_col3 = True
                continue
            elif line == 'col4=endpos':
                in_col3 = False
                continue

            if in_col3:
                try:
                    count, position = map(int, line.split())
                except ValueError:
                    continue  # Skip lines that don't have two integers
                counts[position] = counts.get(position, 0) + count
                total_count += count

    if 51 in counts and counts[51] > total_count / 2:
        return False
    else:
        return True

def format_rec(rec, padlen) -> None:
    '''Process a record and give it a stagger, writing to stdout'''
    # Apply padding
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
@click.argument('info', required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.argument('infosumm', required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.argument('fastq', nargs=2, required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
def stagger_gih(info, infosumm, fastq):
    """
    Given an info file from cutadapt and a processed info_summary file, adds any necessary staggers
    to the barcodes in the input FASTQ file. Writes to stdout. INTENDED FOR INTERNAL USE in harpy preprocess gih.
    """
    stagger = needs_stagger(infosumm)
    t = 0
    discarded = 0
    with (
        open(info, 'r') as f,
        pysam.FastxFile(fastq[0], persist=False) as fq_F,
        pysam.FastxFile(fastq[1], persist=False) as fq_R
    ):
        for cutadapt, fq1, fq2 in zip_longest(f, fq_F, fq_R):
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
                raise RuntimeError(f"Error: Read name mismatch at line {t}. Expected:\n{CA_name}\nbut found:\n{fq1.name}")

            if not stagger:
                sys.stdout.write(str(fq1) + '\n')
                if fq2:
                    sys.stdout.write(str(fq2) + "\n")
                continue

            col3 = int(fields[2])
            plen = 7 if (col3 < 51 or col3 > 58) else (col3 - 51)

            format_rec(fq1, plen)
            sys.stdout.write(str(fq2) + '\n')

    sys.stderr.write(f"Reads without ME sequence that were discarded: {discarded}\n")
