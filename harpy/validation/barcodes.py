"""
Processes and validations relating to identifying barcodes and/or linkd-read types
"""

import pysam
import re
import subprocess
from harpy.common.file_ops import is_gzip, safe_read
from harpy.common.printing import print_error
from harpy.common.progress import  harpy_progresspanel, harpy_progressbar

def validate_barcodefile(infile: str, return_len: bool = False, quiet: int = 0, limit: int = 60, gzip_ok: bool = True, haplotag_only: bool = False, check_dups: bool = True) -> None | int:
    """Does validations to make sure it's one length, within a length limit, one per line, and nucleotides"""
    if is_gzip(infile) and not gzip_ok:
        print_error("incorrect format", f"The input file must be in uncompressed format. Please decompress [blue]{infile}[/] and try again.")
    lengths = set()
    nucleotides = {'A','C','G','T'}
    def validate(line_num, bc_line):
        barcode = bc_line.rstrip()
        if len(barcode.split()) > 1:
            print_error("incorrect format", f"There must be one barcode per line, but multiple entries were detected on [bold]line {line_num}[/] in [blue]{infile}[/]")
        if not set(barcode).issubset(nucleotides) or barcode != barcode.upper():
            print_error("incorrect format", f"Invalid barcode format on [bold]line {line_num }[/]: [yellow]{barcode}[/].\nBarcodes in [blue]{infile}[/] must be captial letters and only contain standard nucleotide characters [green]ATCG[/].")
        return len(barcode)
    progress = harpy_progressbar(quiet)
    with safe_read(infile) as bc_file, harpy_progresspanel(progress, title= "Validating barcodes", quiet=quiet):
        out = subprocess.Popen(['wc', '-l', infile], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
        linenum = int(out.partition(b' ')[0])
        if linenum > 96**4 and haplotag_only:
            print_error("Too many barcodes", f"The maximum number of barcodes possible with haplotagging is [bold]96^4 (84,934,656)[/], but there are [yellow]{linenum}[/] barcodes in [blue]{infile}[/]. Please use fewer barcodes.")
        task_progress = progress.add_task("[dim]Progress", total=linenum)
        # check for duplicates
        if check_dups:
            sort_out = subprocess.Popen(["sort", infile], stdout=subprocess.PIPE)
            dup_out = subprocess.run(["uniq", "-d"], stdin=sort_out.stdout, capture_output=True, text=True)
            if dup_out.stdout:
                print_error(
                    "duplicate barcodes",
                    f"Duplicate barcodes were detected in {infile}, which will result in misleading simulated data.",
                    "Check that you remove duplicate barcodes from your input file.",
                    "Duplicates identified",
                    dup_out.stdout
                )
        for line,bc in enumerate(bc_file, 1):
            length = validate(line, bc)
            if length > limit:
                print_error("barcodes too long", f"Barcodes in [blue]{infile}[/] are [yellow]{length}bp[/] and cannot exceed a length of [bold]{limit}bp[/]. Please use shorter barcodes.")
            lengths.add(length)
            if len(lengths) > 1:
                str_len = ", ".join(str(_length) for _length in lengths)
                print_error("inconsistent length", f"Barcodes in [blue]{infile}[/] must all be a single length, but multiple lengths were detected: [yellow]{str_len}[/]")
            progress.advance(task_progress)
    if not lengths:
        print_error("no barcodes detected", f"No barcodes were found in [blue]{infile}[/]. Please check the input file.")
    if return_len:
        return lengths.pop()

def which_linkedread(fastq: str) -> str:
    """
    Scans the first 100 records of a FASTQ file and tries to determine the barcode technology
    Returns one of: "haplotagging", "stlfr", "tellseq", or "none"
    """
    haplotagging = re.compile(r'\s?BX:Z:(A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2})')
    stlfr = re.compile(r'#([0-9]+_[0-9]+_[0-9]+)(\s|$)')
    tellseq = re.compile(r':([ATCGN]+)(\s|$)')
    recs = 1
    with pysam.FastxFile(fastq, persist=False) as fq:
        for i in fq:
            if recs > 100:
                break
            if i.comment and haplotagging.search(i.comment):
                return "haplotagging"
            elif stlfr.search(i.name):
                return "stlfr"
            elif tellseq.search(i.name):
                return "tellseq"
            recs += 1
    return "none"

def which_linkedread_sam(file_path: str) -> str:
    """
    Scans the first 100 records of a SAM/BAM file and tries to determine the barcode technology
    Returns one of: "haplotagging", "stlfr", "tellseq", or "none"
    """
    recs = 1
    with pysam.AlignmentFile(file_path, require_index=False) as alnfile:
        for record in alnfile.fetch(until_eof = True):
            if recs > 100:
                break
            try:
                bx = record.get_tag("BX")
                if re.search(r"^[ATCGN]+$", bx):
                    return "tellseq"
                elif re.search(r"^\d+_\d+_\d+$", bx):
                    return "stlfr"
                elif re.search(r"^A\d{2}C\d{2}B\d{2}D\d{2}$", bx):
                    return "haplotagging"
                else:
                    recs += 1
                    continue
            except KeyError:
                recs += 1
                continue
    return "none"

