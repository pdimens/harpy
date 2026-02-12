import argparse
import os
import pysam
import shutil
import subprocess
import sys
from itertools import zip_longest

class SimpleCachedFQWriter():
    def __init__(self, prefix: str, cachemax: int = 10000, n_out:int = 2, threads: int = 1):
        """
        A cache for R1 and R2 reads that writes to the stdin of an open gzip subprocesses when the cache exceeds cachemax entries.
        A somewhat more naiive version than the one in Djinn.
        """
        if threads // 2 <= 1:
            compressor = ['pigz', '-p', '1'] if shutil.which("pigz") else ['gzip']
        else:
            compressor = ['pigz', '-p', f"{threads // n_out}"] if shutil.which("pigz") else ['gzip']

        self.R1_out = open(f"{prefix}.R1.fq.gz", "wb")
        self.gz_R1 = subprocess.Popen(compressor, stdout= self.R1_out, stdin=subprocess.PIPE)
        if n_out == 2:
            self.R2_out = open(f"{prefix}.R2.fq.gz", "wb") 
            self.gz_R2 = subprocess.Popen(compressor, stdout= self.R2_out, stdin=subprocess.PIPE)
        else:
            self.R2_out = None
            self.gz_R2 = None

        self.r1_cache: list[bytes] = []
        self.r2_cache: list[bytes] = []
        self.cachemax: int = cachemax

    def queue(self, *records: pysam.FastxRecord|None):
        """add r1 and r2 reads as bytestrings to the cache, write and empty cache when cache exceeds self.max reads"""
        for idx,i in enumerate(records):
            if i:
                j = str(i).encode("utf-8")
                if idx == 0:
                    self.r1_cache.append(j)
                else:
                    self.r2_cache.append(j)

        if len(self.r1_cache) >= self.cachemax or len(self.r2_cache) >= self.cachemax:
            self.write()

    def write(self):
        """writes r1 and r2 cache to self.gq+_R* and empties caches"""
        self.gz_R1.stdin.write(b"\n".join(self.r1_cache))
        self.r1_cache = []
        if self.gz_R2:
            self.gz_R2.stdin.write(b"\n".join(self.r2_cache))
            self.r2_cache = []

    def __enter__(self):
        return self

    def close(self):
        '''flush remaining cache and clean everything up'''
        self.write()

        self.gz_R1.stdin.close()
        self.gz_R1.wait()
        self.R1_out.close()
        if self.gz_R2:
            self.gz_R2.stdin.close()
            self.gz_R2.wait()
            self.R2_out.close()

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.close()


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
    parser.add_argument('prefix', type = str, help = "output filename prefix")   
    parser.add_argument('info', type = str, help = "info file from cutadapt")
    parser.add_argument('fastq', nargs=2, type = str, help = "forwards and reverse fastq files")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if not os.path.isfile(args.info):
        parser.error(f"{args.info} does not exist")
    for i in args.fastq:
        if not os.path.isfile(i):
            parser.error(f"{i} does not exist")

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

    discarded = 0
    try:
        with (
            open(args.info, 'r') as f,
            pysam.FastxFile(args.fastq[0], persist=True) as fq,
            pysam.FastxFile(args.fastq[1], persist=True) as fq2,
            SimpleCachedFQWriter(args.prefix, args.batchsize, n_out = 2, threads = args.threads) as writer
        ):
            t = 0
            for line, entry, entry2 in zip_longest(f, fq, fq2):
                if line and not entry:
                    RuntimeError(f"Error: More reads in {args.info} than in {args.fastq[0]}.")
                if entry and not line:
                    RuntimeError(f"Error: More reads in {args.fastq[0]} than in {args.info}.")
                if not line and not entry and entry2:
                    writer.queue(None, entry2)
                    continue

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
                    entry.sequence = f"{pad[p]}{entry.sequence[1:]}"
                    entry.quality = f"{qpad[p]}{entry.quality[1:]}"
                else:
                    entry.sequence = f"{pad[p]}{entry.sequence}"
                    entry.quality = f"{qpad[p]}{entry.quality}"

                writer.queue(entry, entry2)

    except Exception as e:
        print(f"{e}", file = sys.stderr)

    finally:
        sys.stderr.write(f"Reads discarded without ME sequence: {discarded}\n")
