import argparse
import os
import pysam
import subprocess
import shutil
import sys


def needs_stagger(filename):
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

def process_info_row(line: str):

    fields = line.rstrip().split('\t')
    
    # Extract read ID (first token of first field)
    sp = fields[0].find(' ')
    if sp != -1:
        expID = "@" + fields[0][:sp]
    else:
        expID = "@" + fields[0]
    
    # Compute pad length p
    col2 = int(fields[1])
    col3 = int(fields[2])
    
    if col2 == -1 or col3 < 51 or col3 > 58:
        padlen = 7
    else:
        padlen = col3 - 51

    return expID, padlen

def main():
    parser = argparse.ArgumentParser(
        prog = 'gih_stagger',
        description =
        """
        Determines if fastq files need a stagger for linked-read barcodes and outputs FASTQ records with the appropriate
        stagger if they do. CAUTION: Intended for internal use for GIH-haplotagging chemistry only!
        """,
        usage = "gih_stagger infofile infosummary fastq > output.fq.gz",
        exit_on_error = False
        )

    parser.add_argument('info_file', help = "Info file written by cutadapt --info")
    parser.add_argument('info_summary', help = "Harpy-produced summary file for info_file")
    parser.add_argument('fastq_file', help = "Input FASTQ file")
    parser.add_argument('threads', help = "Number of threads if pigz is on the system")
    parser.add_argument('batchsize', help = "Batch-write this many output fastq records at a time")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    for i in [args.info_file, args.info_summary, args.fastq_file]:
        if not os.path.exists(i):
            parser.error(f"{i} was not found")

    if not needs_stagger(args.info_summary):
        subprocess.run(['cat', args.fastq_file])
        sys.exit(0)

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
    _cmd = ['pigz', '-p', f'{max(args.threads-1, 1)}'] if shutil.which("pigz") else ['gzip']
    batch = []
    
    with (
        open(args.info_file, 'r') as f,
        pysam.FastxFile(args.fastq_file) as fq,
        subprocess.Popen(_cmd, stdin=subprocess.PIPE, stdout=sys.stdout, bufsize=1024*1024) as compressor
    ):
        current_read = 0
        while True:
            current_read += 1
            try:
                line = f.readline()
                expID, padlen = process_info_row(line)
                record = fq.__next__()
            except StopIteration:
                break

            # Extract header token up to first space
            hdr_full = f"@{record.name}"
            sp = hdr_full.find(' ')
            if sp != -1:
                hdr = hdr_full[:sp]
            else:
                hdr = hdr_full
            
            if hdr != expID:
                RuntimeError(f"Error: Read name mismatch at read {current_read}! Expected {expID} but found {hdr}.")
            
            # Build padded FASTQ record
            # Header line (include full comment if present)
            fq_rec = []
            if record.comment:
                fq_rec.append(f"@{record.name}\t{record.comment}")
            else:
                fq_rec.append(f"@{record.name}")
            
            # Sequence line
            if p == 6:
                fq_rec.append(f"{pad[padlen]}{record.sequence[1:]}")
            else:
                fq_rec.append(f"{pad[padlen]}{record.sequence}")
            
            # '+' line
            fq_rec.append("+\n")
            
            # Quality line
            if p == 6:
                fq_rec.append(f"{qpad[padlen]}{record.quality[1:]}")
            else:
                fq_rec.append(f"{qpad[padlen]}{record.quality}")
            
            batch.append("\n".join(fq_rec))
            # Write batch if it reaches batch_size lines
            if len(batch) >= args.batchsize:
                compressor.stdin.write(''.join(batch).encode())
                batch = []
        
        # Write remaining batch
        if batch:
            compressor.stdin.write('\n'.join(batch).encode())
        compressor.stdin.close()
        compressor.wait()    



