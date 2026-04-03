"""
Processes and validations relating to identifying barcodes and/or linked-read types
"""

import subprocess
from harpy.common.file_ops import is_gzip, safe_read
from harpy.common.printing import HarpyPrint

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




