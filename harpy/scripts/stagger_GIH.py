import argparse
import os
import pysam
import shutil
import subprocess
import sys
from itertools import zip_longest

def main():
    parser = argparse.ArgumentParser(
        prog = 'stagger-GIH',
        description =
        """
        Given an info file from cutadapt and a processed info_summary file, adds any necessary staggers
        to the barcodes in the input FASTQ file. Writes to stdout. INTENDED FOR INTERNAL USE in harpy preprocess gih.
        """,
        usage = "stagger-GIH -t 4 -b 20000 INFOFILE FASTQ > OUT_FASTQ",
        exit_on_error = False
        )

    parser.add_argument('-t', '--threads', type = int, default = 4, help = "Number of threads for pigz (if available)")
    parser.add_argument('-b', '--batchsize', type = int, default = 10000, help = "Batch write this many FASTQ records at a time")
    parser.add_argument('info', type = str, help = "info file from cutadapt")
    parser.add_argument('fastq', type = str, help = "fastq file")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if not os.path.isfile(args.info):
        parser.error(f"{args.info} does not exist")
    if not os.path.isfile(args.fastq):
        parser.error(f"{args.fastq} does not exist")

    # Precompute pad sequences for p=0..7
    pad = {
        0: "TTTTTTT",
        1: "CCCCCC",
        2: "GGGGG",
        3: "AAAA",
        4: "TTT",
        5: "CC",
        6: "GG",
        7: ""
    }
    
    # Precompute quality pads (7-p) × 'I'
    qpad = {}
    for i in range(8):
        qpad[i] = "I" * (7 - i)

    # Setup compression subprocess
    _cmd = ['pigz', '-c', '-p', f'{max(args.threads-1, 1)}'] if shutil.which("pigz") else ['gzip', '-c']
    compressor = subprocess.Popen(
        _cmd,
        stdin=subprocess.PIPE,
        stdout=sys.stdout,
        stderr=sys.stderr,
        bufsize=1024*1024
    )

    batch = []
    discarded = 0
    try:
        with open(args.info, 'r') as f, pysam.FastxFile(args.fastq, persist=False) as fq:
            t = 0
            for line, entry in zip_longest(f, fq):
                if line and not entry:
                    RuntimeError(f"Error: More reads in {args.info} than in {args.fastq}.")
                if entry and not line:
                    RuntimeError(f"Error: More reads in {args.fastq} than in {args.info}.")

                t += 1
                line = line.rstrip('\n')
                fields = line.split('\t')

                # Extract read ID (first token of first field)
                n = fields[0].split(' ')[0]

                # Compute pad length p
                col2 = int(fields[1])
                if col2 == -1:
                    discarded += 1
                    continue
                col3 = int(fields[2])
                p = 7 if (col3 < 51 or col3 > 58) else (col3-51)

                _name = f"{entry.name}"
                if _name != n:
                    RuntimeError(f"Error: Read name mismatch at line {t}! Expected {n} in fastq but found {_name}.")

                # Build padded FASTQ entry
                if entry.comment:
                    _name += f"\t{entry.comment}"

                if p == 6:
                    _seq = f"{pad[p]}{entry.sequence[1:]}"
                    _qual = f"{qpad[p]}{entry.quality[1:]}"
                else:
                    _seq = f"{pad[p]}{entry.sequence}"
                    _qual = f"{qpad[p]}{entry.quality}"

                batch.append(
                    f"@{_name}\n"
                    f"{_seq}\n"
                    f"+\n"
                    f"{_qual}\n"
                )
                # Write batch if it reaches batch_size lines
                if len(batch) >= args.batchsize:
                    compressor.stdin.write(''.join(batch).encode())
                    batch = []
        
        # Write remaining batch
        if batch:
            compressor.stdin.write(''.join(batch).encode())
    
    except Exception as e:
        print(f"{e}", file = sys.stderr)
        compressor.kill()
    
    finally:
        sys.stderr.write(f"Reads discarded without ME sequence: {discarded}\n")
        compressor.stdin.close()
        compressor.wait()    
