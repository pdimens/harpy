"""
Processes and validations relating to identifying barcodes and/or linked-read types
"""

import pysam
import re
import subprocess
import shutil
from harpy.common.file_ops import is_gzip, safe_read
from harpy.common.printing import HarpyPrint

HAPLOTAGGING_RX = re.compile(r'\s?BX:Z:(A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2})')
HAPLOTAGGING_RX_SAM = re.compile(r"^A\d{2}C\d{2}B\d{2}D\d{2}$")
STLFR_RX = re.compile(r'#([0-9]+_[0-9]+_[0-9]+)(\s|$)')
STLFR_RX_SAM = re.compile(r"^\d+_\d+_\d+$")
TELLSEQ_RX = re.compile(r':([ATCGN]+)(\s|$)')
TELLSEQ_RX_SAM = re.compile(r"^[ATCGN]+$")

def validate_barcodefile(infile: str, return_len: bool = False, quiet: int = 0, limit: int = 60, gzip_ok: bool = True, haplotag_only: bool = False, check_dups: bool = True) -> None | int:
    """Does validations to make sure it's one length, within a length limit, one per line, and nucleotides"""
    hp = HarpyPrint(quiet)
    if is_gzip(infile) and not gzip_ok:
        hp.error("incorrect format", f"The input file must be in uncompressed format. Please decompress [blue]{infile}[/] and try again.")
    lengths = set()
    nucleotides = {'A','C','G','T'}
    def validate(line_num, bc_line):
        barcode = bc_line.rstrip()
        if len(barcode.split()) > 1:
            hp.error("incorrect format", f"There must be one barcode per line, but multiple entries were detected on [bold]line {line_num}[/] in [blue]{infile}[/]")
        if not set(barcode).issubset(nucleotides) or barcode != barcode.upper():
            hp.error("incorrect format", f"Invalid barcode format on [bold]line {line_num }[/]: [yellow]{barcode}[/].\nBarcodes in [blue]{infile}[/] must be captial letters and only contain standard nucleotide characters [green]ATCG[/].")
        return len(barcode)
    progress = hp.progressbar()
    with safe_read(infile) as bc_file, hp.progresspanel(progress, title= "Validating barcodes"):
        out = subprocess.Popen(['wc', '-l', infile], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
        linenum = int(out.partition(b' ')[0])
        if linenum > 96**4 and haplotag_only:
            hp.error("Too many barcodes", f"The maximum number of barcodes possible with haplotagging is [bold]96^4 (84,934,656)[/], but there are [yellow]{linenum}[/] barcodes in [blue]{infile}[/]. Please use fewer barcodes.")
        task_progress = progress.add_task("[dim]Progress", total=linenum)
        # check for duplicates
        if check_dups:
            sort_out = subprocess.Popen(["sort", infile], stdout=subprocess.PIPE)
            dup_out = subprocess.run(["uniq", "-d"], stdin=sort_out.stdout, capture_output=True, text=True)
            if dup_out.stdout:
                hp.error(
                    "duplicate barcodes",
                    f"Duplicate barcodes were detected in {infile}, which will result in misleading simulated data.",
                    "Check that you remove duplicate barcodes from your input file.",
                    "Duplicates identified",
                    dup_out.stdout
                )
        for line,bc in enumerate(bc_file, 1):
            length = validate(line, bc)
            if length > limit:
                hp.error("barcodes too long", f"Barcodes in [blue]{infile}[/] are [yellow]{length}bp[/] and cannot exceed a length of [bold]{limit}bp[/]. Please use shorter barcodes.")
            lengths.add(length)
            if len(lengths) > 1:
                str_len = ", ".join(str(_length) for _length in lengths)
                hp.error("inconsistent length", f"Barcodes in [blue]{infile}[/] must all be a single length, but multiple lengths were detected: [yellow]{str_len}[/]")
            progress.advance(task_progress)
    if not lengths:
        hp.error("no barcodes detected", f"No barcodes were found in [blue]{infile}[/]. Please check the input file.")
    if return_len:
        return lengths.pop()


def which_linkedread(fastq: str) -> str:
    """
    Scans the first 100 records of a FASTQ file and tries to determine the barcode technology
    Returns one of: "haplotagging", "stlfr", "tellseq", or "none"
    """
    with pysam.FastxFile(fastq, persist=False) as fq:
        for i,record in enumerate(fq, 1):
            if i > 100:
                break
            if record.comment and HAPLOTAGGING_RX.search(record.comment):
                return "haplotagging"
            if STLFR_RX.search(record.name):
                return "stlfr"
            if TELLSEQ_RX.search(record.name):
                return "tellseq"
    return "none"


def which_linkedread_sam(file_path: str) -> str:
    """
    Scans the first 100 records of a SAM/BAM file and tries to determine the barcode technology
    Returns one of: "haplotagging", "stlfr", "tellseq", or "none"
    """
    with pysam.AlignmentFile(file_path, require_index=False) as alnfile:
        for i, record in enumerate(alnfile.fetch(until_eof = True), 1):
            if i > 100:
                break
            if not record.has_tag("BX"):
                continue
            bx = record.get_tag("BX")
            if TELLSEQ_RX_SAM.search(bx):
                return "tellseq"
            if STLFR_RX_SAM.search(bx):
                return "stlfr"
            if HAPLOTAGGING_RX_SAM.search(bx):
                return "haplotagging"
    return "none"

#def is_staggered(infile: str) -> bool:
#    count_51 = 0
#    total_count = 0
#    in_col3 = False
#
#    with open(infile, 'r') as file:
#        for line in file:
#            line = line.strip()
#            
#            if line == 'col3=startpost':
#                in_col3 = True
#                continue
#            elif line == 'col4=endpos':
#                break  # No need to read the rest of the file
#            
#            if in_col3:
#                try:
#                    count, position = line.split()
#                    count = int(count)
#                    if position == '51':  # String comparison before conversion
#                        count_51 += count
#                    total_count += count
#                except ValueError:
#                    continue
#
#    return count_51 <= total_count / 2
#
#
#def parse_infofile(info_file: str):
#    """Parse sample_info.txt and return lists of expected IDs and pad lengths."""
#    exp_ids = []
#    pad_lens = []
#    
#    with open(info_file) as f:
#        for line in f:
#            fields = line.strip().split('\t')
#            rid = fields[0].split()[0]
#            
#            try:
#                col2 = int(fields[1])
#                col3 = int(fields[2])
#                p = 7 if (col2 == -1 or col3 < 51 or col3 > 58) else col3 - 51
#            except (IndexError, ValueError):
#                p = 7
#            
#            exp_ids.append(rid)
#            pad_lens.append(p)
#    
#    return exp_ids, pad_lens

def stagger_info(info_file: str):
    """
    Parse sample_info.txt. Returns (exp_ids: list[str], pad_lens: list[int]) if staggered, None otherwise
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
