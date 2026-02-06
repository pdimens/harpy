"""
Functions for GIH-haplotagging preprocessing
"""

import pysam
import subprocess
import shutil

def stagger_info(info_file: str):
    """
    Parse `sample_info.txt` from `cutadapt`. Returns (exp_ids: `list[str]`, pad_lens: `list[int]`) if staggered, `None` otherwise.
    """
    exp_ids = []
    pad_lens = []
    count_51 = 0
    total_count = 0
    in_col3 = False
    
    with open(info_file) as f:
        for line in f:
            line = line.strip()

            # Check for section markers
            if line == 'col3=startpost':
                in_col3 = True
                continue
            elif line == 'col4=endpos':
                in_col3 = False
                continue

            # Process col3 section for stagger detection
            if in_col3:
                try:
                    count, position = line.split()
                    count = int(count)
                    if position == '51':
                        count_51 += count
                    total_count += count
                except ValueError:
                    pass
                continue

            # Parse main data (before col3 section or after col4)
            if not in_col3 and line and not line.startswith('col'):
                try:
                    fields = line.split('\t')
                    rid = fields[0].split()[0]

                    col2 = int(fields[1])
                    col3 = int(fields[2])
                    p = 7 if (col2 == -1 or col3 < 51 or col3 > 58) else col3 - 51

                    exp_ids.append(rid)
                    pad_lens.append(p)
                except (IndexError, ValueError):
                    # Skip malformed lines
                    pass

    # Check if staggered
    is_staggered = count_51 <= total_count / 2

    if is_staggered:
        return exp_ids, pad_lens
    else:
        return None

def process_stagger(fastq_file :str, output_file: str, exp_ids: list[str], pad_lens: list[int], threads: int, batchsize: int = 10000):
    """Process FASTQ stagger using pysam with batched compression."""
    # Precompute pads
    pad = ["TTTTTTT", "CCCCCC", "GGGGG", "AAAA", "TTT", "CC", "GG", ""]
    qpad = ["IIIIIII", "IIIIII", "IIIII", "IIII", "III", "II", "I", ""]

    BATCH_SIZE = batchsize
    _cmd = ['pigz', '-c', '-p', f'{max(threads-1, 1)}'] if shutil.which("pigz") else ['gzip', '-c']
    compressor = subprocess.Popen(
        _cmd,
        stdin=subprocess.PIPE,
        stdout=open(output_file, 'wb'),
        bufsize=1024*1024
    )
    batch = []

    try:
        with pysam.FastxFile(fastq_file, persist=False) as fq:
            for read_idx, read in enumerate(fq):
                read_name = read.name.split()[0]
                if read_name != exp_ids[read_idx]:
                    compressor.stdin.close()
                    compressor.wait()
                    raise ValueError(
                        f"Read name mismatch at read {read_idx+1}. "
                        f"Expected {exp_ids[read_idx]} but found {read_name}\n Terminating.\n"
                    )
                p = pad_lens[read_idx]

                # Apply padding
                if p == 6:
                    padded_seq = pad[p] + read.sequence[1:]
                else:
                    padded_seq = pad[p] + read.sequence
                
                padded_qual = qpad[p] + read.quality

                # Don't store 'read' object, as it becomes invalid after next iteration
                record = f"@{read.name}\n{padded_seq}\n+\n{padded_qual}\n"
                batch.append(record)

                if len(batch) >= BATCH_SIZE:
                    compressor.stdin.write(''.join(batch).encode('utf-8'))
                    batch = []

        if batch:
            compressor.stdin.write(''.join(batch).encode('utf-8'))

    finally:
        compressor.stdin.close()
        compressor.wait()
