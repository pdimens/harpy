"""
Functions for GIH-haplotagging preprocessing
"""

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


def stagger_info(info_file: str):
    """
    Parse `sample_info.txt` from `cutadapt`. Returns (exp_ids: `list[str]`, pad_lens: `list[int]`) if staggered, `None` otherwise.
    """
    exp_ids = []
    pad_lens = []
    count_51 = 0
    total_count = 0
    in_col3 = False
    
    with open(info_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Check for section markers
            if line == 'col3=startpos':
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

#!/usr/bin/env python3


def padUMI(info_file, fastq_file, fastq_out, logfile, threads, batchsize):
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
    
    # Read sample_info.txt
    expID = {}
    plen = {}
    ninfo = 0
    
    with open(info_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            fields = line.split('\t')
            
            # Extract read ID (first token of first field)
            sp = fields[0].find(' ')
            if sp != -1:
                rid = fields[0][:sp]
            else:
                rid = fields[0]
            
            # Compute pad length p
            col2 = int(fields[1])
            col3 = int(fields[2])
            
            if col2 == -1 or col3 < 51 or col3 > 58:
                p = 7
            else:
                p = col3 - 51
            
            ninfo += 1
            expID[ninfo] = "@" + rid
            plen[ninfo] = p
    
    # Setup compression subprocess
    _cmd = ['pigz', '-c', '-p', f'{max(threads-1, 1)}'] if shutil.which("pigz") else ['gzip', '-c']
    compressor = subprocess.Popen(
        _cmd,
        stdin=subprocess.PIPE,
        stdout=open(fastq_out, 'wb'),
        bufsize=1024*1024
    )
    
    batch = []

    try:
        # Process FASTQ file with pysam
        t = 0  # read counter
        
        with pysam.FastxFile(fastq_file) as fh:
            for entry in fh:
                t += 1
                
                if t > ninfo:
                    RuntimeError(f"Error: More reads in {fastq_file} than in {info_file}. Read {t}")
                
                p = plen[t]
                n = expID[t]
                
                # Extract header token up to first space
                hdr_full = f"@{entry.name}"
                sp = hdr_full.find(' ')
                if sp != -1:
                    hdr = hdr_full[:sp]
                else:
                    hdr = hdr_full
                
                if hdr != n:
                    RuntimeError(f"Error: Read name mismatch at read {t}! Expected {n} but found R1={hdr}.")
                
                # Build padded FASTQ entry
                # Header line (include full comment if present)
                fq_rec = []
                if entry.comment:
                    fq_rec.append(f"@{entry.name}\t{entry.comment}")
                else:
                    fq_rec.append(f"@{entry.name}")
                
                # Sequence line
                if p == 6:
                    fq_rec.append(f"{pad[p]}{entry.sequence[1:]}")
                else:
                    fq_rec.append(f"{pad[p]}{entry.sequence}")
                
                # '+' line
                fq_rec.append("+\n")
                
                # Quality line
                if p == 6:
                    fq_rec.append(f"{qpad[p]}{entry.quality[1:]}")
                else:
                    fq_rec.append(f"{qpad[p]}{entry.quality}")
                
                batch.append("\n".join(fq_rec))
                # Write batch if it reaches batch_size lines
                if len(batch) >= batchsize:
                    compressor.stdin.write(''.join(batch).encode())
                    batch = []
        
        # Write remaining batch
        if batch:
            compressor.stdin.write('\n'.join(batch).encode())
        
        # Check if we processed all expected reads
        if t < ninfo:
            RuntimeError(f"Error: Fewer reads in {t} than in {ninfo}.")
    
    except Exception as e:
        print(f"{e}")
        compressor.kill()
    
    finally:
        compressor.stdin.close()
        compressor.wait()    





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
