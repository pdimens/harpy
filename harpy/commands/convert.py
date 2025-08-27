#! /usr/bin/env python3

from concurrent.futures import ThreadPoolExecutor
from itertools import product, zip_longest
import os
from random import sample
import re
import subprocess
import sys
import rich_click as click
import pysam
from harpy.common.file_ops import safe_read
from harpy.common.cli_filetypes import SAMfile, FASTQfile
from harpy.common.convert import FQRecord, compress_fq, INVALID_10x, INVALID_HAPLOTAGGING, INVALID_STLFR, INVALID_TELLSEQ
from harpy.common.printing import print_error
from harpy.common.progress import harpy_pulsebar
from harpy.validation.barcodes import validate_barcodefile, which_linkedread

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def convert():
    """
    Convert between linked-read formats and barcode styles
    """

module_docstring = {
    "harpy convert": [
        {
            "name": "Commands",
            "commands": ["fastq", "standardize-bam", "standardize-fastq"],
            "panel_styles": {"border_style" : "blue"}
        }
    ]
}

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/convert")
@click.option('-o','--output', type = str, metavar= "PREFIX", help='file prefix for output fastq files', required=True)
@click.option('-b','--barcodes', type = click.Path(exists=True, readable=True, dir_okay=False), help='barcodes file [from 10x only]', required=False)
@click.option('--quiet', show_default = True, default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` `1` (all) or `2` (no) output')
@click.argument('target', metavar = "TARGET", type = click.Choice(["10x", "haplotagging", "stlfr", "tellseq"], case_sensitive=False), nargs = 1)
@click.argument('fq1', metavar="R1_FASTQ", type = FASTQfile(single=True), required = True, nargs=1)
@click.argument('fq2', metavar="R2_FASTQ", type = FASTQfile(single=True), required=True, nargs= 1)
def fastq(target,fq1,fq2,output,barcodes, quiet):
    """
    Convert between linked-read FASTQ formats
    
    Autodetects the input data format and takes the positional argument `TARGET` specifying the target data format.
    10X input data requires a `--barcodes` file containing one nucleotide barcode per line to
    determine which barcodes are valid/invalid. In all cases, a file will be created with
    the barcode conversion map. Requires 2 threads.
    
    | from/to      | barcode format                                     | example                     |
    |:-------------|:---------------------------------------------------|:----------------------------|
    | 10x          | the first N base pairs of R1, given `--barcodes`   |                             |
    | haplotagging | a `BX:Z:ACBD` SAM tag in the sequence header       | `@SEQID BX:Z:A01C93B56D11`  |
    | stlfr        | `#1_2_3` format appended to the sequence ID        | `@SEQID#1_2_3`              |
    | tellseq      | `:ATCG` format appended to the sequence ID         | `@SEQID:GGCAAATATCGAGAAGTC` |
    """
    from_ = which_linkedread(fq1)
    if from_ == "none" and not barcodes:
        print_error("missing barcodes file", "The input data was inferred to be 10X format because neither haplotagging, stlfr, nor tellseq barcodes were detected in their appropriate locations in the first 100 fastq records of the forward reads. A [green]--barcodes[/] file must be provided if the input data is 10x. If the input data is not 10X, then it is formatted incorrectly for whatever technology it was generated with.")
    if from_ != "none" and barcodes:
        print_error("unsupported process", f"The input data was inferred to be {from_} format, but a [green]--barcodes[/] file was provided, suggesting it's 10X format data. If the input data is not {from_}, then it is formatted incorrectly for whatever technology it was generated with.")
    if from_ == target:
        print_error("identical conversion target", f"The input file was inferred to be[green]{from_}[/], which is identical to the conversion target [green]{target}[/]. The formats must be different from each other. If the input data is not {from_}, then it is formatted incorrectly for whatever technology it was generated with.")

    # just make sure it's all lowercase
    to_ = target.lower()
    # for barcodes, use sample() so the barcodes don't all start with AAAAAAAAAAAAA (or 1)
    # it's not functionally important, but it does make the barcodes *look* more distinct
    if to_ == "10x":
        bc_generator = product(*[sample("ATCG", 4) for i in range(16)])
        invalid = INVALID_10x
        def format_bc(bc):
            return "".join(bc)
    elif to_ == "tellseq":
        bc_generator = product(*[sample("ATCG", 4) for i in range(18)])
        invalid = INVALID_TELLSEQ
        def format_bc(bc):
            return "".join(bc)
    elif to_ == "stlfr":
        bc_generator = product(*[sample(range(1,1538), 1537) for i in range(3)])
        invalid = INVALID_STLFR
        def format_bc(bc):
            return "_".join(str(i) for i in bc)
    elif to_ == "haplotagging":
        bc_generator = product(
            ["A" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["C" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["B" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["D" + str(i).zfill(2) for i in sample(range(1,97), 96)]
        )
        invalid = INVALID_HAPLOTAGGING
        def format_bc(bc):
            return "".join(bc)

    if from_ == "10x":
        bc_len = validate_barcodefile(barcodes, return_len=True, quiet=quiet, check_dups=False)
        with safe_read(barcodes) as b:
            barcodelist = set(b.readlines())
    else:
        bc_len = 0
    bc_inventory = {}
    # create the output directory in case it doesn't exist
    if os.path.dirname(output):
        os.makedirs(os.path.dirname(output), exist_ok=True)

    with (
        pysam.FastxFile(fq1, persist=False) as R1,
        pysam.FastxFile(fq2, persist=False) as R2,
        open(f"{output}.R1.fq", "w") as R1_out,
        open(f"{output}.R2.fq", "w") as R2_out,
        open(f"{output}.bc", "w") as bc_out,
        harpy_pulsebar(quiet) as progress,
    ):
        progress.add_task(f"[blue]{from_}[/] -> [magenta]{to_}[/]", total = None)
        for r1,r2 in zip_longest(R1,R2):
            if r1:
                _r1 = FQRecord(r1, True, from_, bc_len)
                # 10x validation
                if from_ == "10x":
                    _r1.valid = _r1.barcode in barcodelist
                if _r1.barcode not in bc_inventory:
                    # if it's just tellseq<->10x, keep the existing nucleotide barcode
                    if (from_ in ["tellseq", "10x"] and to_ in ["tellseq", "10x"]):
                        bc_inventory[_r1.barcode] = _r1.barcode
                    else:
                        if _r1.valid:
                            try:
                                bc_inventory[_r1.barcode] = format_bc(next(bc_generator))
                            except StopIteration:
                                print_error("too many barcodes", f"There are more {from_} barcodes in the input data than it is possible to generate {to_} barcodes from.")
                        else:
                            bc_inventory[_r1.barcode] = invalid
                    bc_out.write(f"{_r1.barcode}\t{bc_inventory[_r1.barcode]}\n")
                converted_bc = bc_inventory[_r1.barcode]
                R1_out.write(str(_r1.convert(to_, converted_bc)))
            if r2:
                if r1 and from_ == "10x":
                    _bc = _r1.barcode
                elif not r1 and from_ == "10x":
                    _bc = INVALID_10x
                else:
                    _bc = from_
                # if input format is 10x, copy the barcode to R2
                _r2 = FQRecord(r2, False, _bc, bc_len)
                # check the inventory for existing barcode match
                if _r2.barcode not in bc_inventory:
                    # if it's just tellseq<->10x, keep the existing nucleotide barcode
                    if (from_ in ["tellseq", "10x"] and to_ in ["tellseq", "10x"]):
                        bc_inventory[_r2.barcode] = _r2.barcode
                    else:
                        if _r2.valid:
                            try:
                                bc_inventory[_r2.barcode] = format_bc(next(bc_generator))
                            except StopIteration:
                                print_error("too many barcodes", f"There are more {from_} barcodes in the input data than it is possible to generate {to_} barcodes from.")
                        else:
                            bc_inventory[_r2.barcode] = invalid
                    bc_out.write(f"{_r2.barcode}\t{bc_inventory[_r2.barcode]}\n")
                converted_bc = bc_inventory[_r2.barcode]
                R2_out.write(str(_r2.convert(to_, converted_bc)))
    
    # bgzip compress the output, one file per thread
    with ThreadPoolExecutor(max_workers=2) as executor:
        executor.submit(compress_fq, f"{output}.R1.fq")
        executor.submit(compress_fq, f"{output}.R2.fq")

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/convert")
@click.option('-s', '--style', type = click.Choice(["haplotagging", "stlfr", "tellseq", "10x"], case_sensitive=False), help = 'Change the barcode style')
@click.option('--quiet', show_default = True, default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` `1` (all) or `2` (no) output')
@click.argument('sam', metavar="BAM", type = SAMfile(single=True), required = True, nargs=1)
def standardize_bam(sam, style, quiet):
    """
    Move barcode to BX:Z/VX:i tags in alignments

    If barcodes are present in the sequence name (stlfr, tellseq), this method moves the barcode to the `BX:Z`
    tag of the alignment, maintaining the same barcode style by default. If moved to or already in a `BX:Z` tag,
    will then write a complementary `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid).
    Use `--style` to also convert the barcode to a different style (`haplotagging`, `stlfr`, `tellseq`, `10X`),
    which also writes a `conversion.bc` file to the working directory mapping the barcode conversions. Writes to `stdout`.

    | Option         | Style                                        |
    |:---------------|:---------------------------------------------|
    | `haplotagging` | AxxCxxBxxDxx                                 |
    | `stlfr`        | 1_2_3                                        |
    | `tellseq`      | 18-base nucleotide (e.g. AGCCATGTACGTATGGTA) |
    | `10X`          | 16-base nucleotide (e.g. GGCTGAACACGTGCAG)   |
    """
    logtext = f"Standardizing [dim][-> [magenta]{style.lower()}[/]][/]" if style else "Standardizing"
    convert = None
    if style:
        bc_out = open("conversion.bc", "w")
        convert = style.lower()
        if convert == "tellseq":
            bc_generator = product(*[sample("ATCG", 4) for i in range(18)])
            invalid = INVALID_TELLSEQ
            def format_bc(bc):
                return "".join(bc)
        if convert == "10x":
            bc_generator = product(*[sample("ATCG", 4) for i in range(16)])
            invalid = INVALID_10x
            def format_bc(bc):
                return "".join(bc)
        elif convert == "stlfr":
            bc_generator = product(*[sample(range(1,1538), 1537) for i in range(3)])
            invalid = INVALID_STLFR
            def format_bc(bc):
                return "_".join(str(i) for i in bc)
        elif convert == "haplotagging":
            bc_generator = product(
                ["A" + str(i).zfill(2) for i in sample(range(1,97), 96)],
                ["C" + str(i).zfill(2) for i in sample(range(1,97), 96)],
                ["B" + str(i).zfill(2) for i in sample(range(1,97), 96)],
                ["D" + str(i).zfill(2) for i in sample(range(1,97), 96)]
            )
            invalid = INVALID_HAPLOTAGGING
            def format_bc(bc):
                return "".join(bc)
        bc_inventory = {}

    with (
        pysam.AlignmentFile(sam, require_index=False) as samfile, 
        pysam.AlignmentFile(sys.stdout, "wb", template=samfile) as outfile,
        harpy_pulsebar(quiet, True) as progress,
    ):
        progress.add_task(logtext, total = None)
        for record in samfile.fetch(until_eof=True):
            if record.has_tag("BX") and record.has_tag("VX"):
                print_error("BX/VX tags present", f"The BX:Z and VX:i tags are already present in {os.path.basename(sam)} and does not need to be standardized.")
            if record.has_tag("BX"):
                bx_sanitized = record.get_tag("BX")
                if "0" in bx_sanitized.split("_") or re.search(r"(?:N|[ABCD]00)", bx_sanitized):
                    record.set_tag("VX", 0, "i")
                    vx = 0
                else:
                    record.set_tag("VX", 1, "i")
                    vx = 1
            else:
                # matches either tellseq or stlfr   
                bx = re.search(r"(?:\:[ATCGN]+$|#\d+_\d+_\d+$)", record.query_name)
                if bx:
                    # the 1:0 ignores the first character, which will either be : or #
                    bx_sanitized = bx[0][1:]
                    record.query_name = record.query_name.remove_suffix(bx[0])
                    if "0" in bx_sanitized.split("_") or "N" in bx_sanitized:
                        record.set_tag("VX", 0, "i")
                        vx = 0
                    else:
                        record.set_tag("VX", 1, "i")
                        vx = 1
                    record.set_tag("BX", bx_sanitized, "Z")
                else:
                    outfile.write(record)
            if convert:
                if bx_sanitized not in bc_inventory:
                    if bool(vx):
                        try:
                            bc_inventory[bx_sanitized] = format_bc(next(bc_generator))
                        except StopIteration:
                            print_error("too many barcodes", f"There are more barcodes in the input data than it is possible to generate {convert} barcodes from.")
                    else:
                        bc_inventory[bx_sanitized] = invalid
                    bc_out.write(f"{bx_sanitized}\t{bc_inventory[bx_sanitized]}\n")
                    record.set_tag("BX", bc_inventory[bx_sanitized], "Z")
            outfile.write(record)
    if convert:
        bc_out.close()

@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/convert")
@click.option('--quiet', show_default = True, default = 0, type = click.IntRange(0,2,clamp=True), help = '`0` `1` (all) or `2` (no) output')
@click.option('-s', '--style', type = click.Choice(["haplotagging", "stlfr", "tellseq", "10x"], case_sensitive=False), help = 'Change the barcode style')
@click.argument('prefix', metavar="output_prefix", type = str, required = True, nargs=1)
@click.argument('r1_fastq', metavar="R1_fastq", type = FASTQfile(single=True), required = True, nargs=1)
@click.argument('r2_fastq', metavar="R2_fastq", type = FASTQfile(single=True), required = True, nargs=1)
def standardize_fastq(prefix, r1_fastq, r2_fastq, style, quiet):
    """
    Move barcodes to BX:Z/VX:i tags in sequence headers

    This conversion moves the barcode to the `BX:Z` tag in fastq records, maintaining the same barcode type by default (auto-detected).
    See `harpy convert fastq` for the location and format expectations for different linked-read technologies.
    Also writes a `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid). Use `harpy convert fastq`
    if your data is in 10X format, as this command will not work on 10X format (i.e. barcode is the first 16 bases of read 1).
    Use `--style` to also convert the barcode to a different style (`haplotagging`, `stlfr`, `tellseq`, `10X`).

    | Option         | Style                                        |
    |:---------------|:---------------------------------------------|
    | `haplotagging` | AxxCxxBxxDxx                                 |
    | `stlfr`        | 1_2_3                                        |
    | `tellseq`      | 18-base nucleotide (e.g. AGCCATGTACGTATGGTA) |
    | `10X`          | 16-base nucleotide (e.g. GGCTGAACACGTGCAG)   |
    """
    logtext = f"Standardizing [dim][-> [magenta]{style.lower()}[/]][/]" if style else "Standardizing"
    convert = None
    BC_TYPE = which_linkedread(r1_fastq)
    if not BC_TYPE:
        print_error("undertermined file type", f"Unable to determine the linked-read barcode type after scanning the first 100 records of {os.path.basename(r1_fastq)}. Please make sure the format is one of haplotagging, stlfr, or tellseq. 10X-style with the barcode as the first 16 nucleotides of read 1 is not supported here.")
    
    # create the output directory in case it doesn't exist
    if os.path.dirname(prefix):
        os.makedirs(os.path.dirname(prefix), exist_ok=True)

    if style:
        bc_out = open(f"{prefix}.bc", "w")
        convert = style.lower()
        if convert == "tellseq":
            bc_generator = product(*[sample("ATCG", 4) for i in range(18)])
            invalid = INVALID_TELLSEQ
            def format_bc(bc):
                return "".join(bc)
        if convert == "10x":
            bc_generator = product(*[sample("ATCG", 4) for i in range(16)])
            invalid = INVALID_10x
            def format_bc(bc):
                return "".join(bc)
        elif convert == "stlfr":
            bc_generator = product(*[sample(range(1,1538), 1537) for i in range(3)])
            invalid = INVALID_STLFR
            def format_bc(bc):
                return "_".join(str(i) for i in bc)
        elif convert == "haplotagging":
            bc_generator = product(
                ["A" + str(i).zfill(2) for i in sample(range(1,97), 96)],
                ["C" + str(i).zfill(2) for i in sample(range(1,97), 96)],
                ["B" + str(i).zfill(2) for i in sample(range(1,97), 96)],
                ["D" + str(i).zfill(2) for i in sample(range(1,97), 96)]
            )
            invalid = INVALID_HAPLOTAGGING
            def format_bc(bc):
                return "".join(bc)
        bc_inventory = {}

    with (
        pysam.FastxFile(r1_fastq, persist=False) as R1,
        pysam.FastxFile(r2_fastq, persist=False) as R2,
        open(f"{prefix}.R1.fq", "w") as R1_out,
        open(f"{prefix}.R2.fq", "w") as R2_out,
        harpy_pulsebar(quiet, True) as progress,
    ):
        progress.add_task(logtext, total = None)
        for r1,r2 in zip_longest(R1,R2):
            if r1:
                _r1 = FQRecord(r1, True, BC_TYPE, 0)
                if convert:
                    if _r1.barcode not in bc_inventory:
                        if _r1.valid:
                            try:
                                bc_inventory[_r1.barcode] = format_bc(next(bc_generator))
                            except StopIteration:
                                print_error("too many barcodes", f"There are more {BC_TYPE} barcodes in the input data than it is possible to generate {convert} barcodes from.")
                        else:
                            bc_inventory[_r1.barcode] = invalid
                        # write the barcode to file
                        bc_out.write(f"{_r1.barcode}\t{bc_inventory[_r1.barcode]}\n")
                    # overwrite the record's barcode
                    _r1.barcode = bc_inventory[_r1.barcode]
                R1_out.write(str(_r1.convert("standard", _r1.barcode)))
            if r2:
                _r2 = FQRecord(r2, False, BC_TYPE, 0)
                if convert:
                    if _r2.barcode not in bc_inventory:
                        if _r2.valid:
                            try:
                                bc_inventory[_r2.barcode] = format_bc(next(bc_generator))
                            except StopIteration:
                                print_error("too many barcodes", f"There are more {BC_TYPE} barcodes in the input data than it is possible to generate {convert} barcodes from.")
                        else:
                            bc_inventory[_r2.barcode] = invalid
                        # write the barcode to file
                        bc_out.write(f"{_r2.barcode}\t{bc_inventory[_r2.barcode]}\n")
                    # overwrite the record's barcode
                    _r2.barcode = bc_inventory[_r2.barcode]
                R2_out.write(str(_r2.convert("standard", _r2.barcode)))
    if convert:
        bc_out.close()
        # bgzip compress the output, one file per thread
    with ThreadPoolExecutor(max_workers=2) as executor:
        executor.submit(compress_fq, f"{prefix}.R1.fq")
        executor.submit(compress_fq, f"{prefix}.R2.fq")

@click.command(hidden = True, no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/ncbi")
@click.option('-m', '--barcode-map',  is_flag = True, default = False, help = 'Write a map of the barcode-to-nucleotide conversion')
@click.option('-p', '--prefix', required=True, type = str, help = "Output file name prefix")
@click.option('-i', '--preserve-invalid',  is_flag = True, default = False, help = 'Retain the uniqueness of invalid barcodes')
@click.option('-s', '--scan', show_default = True, default = 100, type = click.IntRange(min=1), help = 'Number of reads to scan to identify barcode location and format')
@click.argument('r1_fq', required=True, type=FASTQfile(single=True), nargs=1)
@click.argument('r2_fq', required=True, type=FASTQfile(single=True), nargs=1)
def ncbi(prefix, r1_fq, r2_fq, scan, preserve_invalid, barcode_map):
    """
    Convert FASTQ files for NCBI submission

    Converts the linked-read barcodes into nucleotide format (if necessary) and adds it to the beginning
    of the sequence, retaining the linked-read barcodes after NCBI/SRA reformats the FASTQ after submission.
    The barcode will be stored as the first 18 bases in both R1 and R2 reads, followed by a 5bp spacer of "NNNNN",
    then the actual sequence. Requires a `--prefix` to name the output files. Use `--barcode-map`/`-m` to write a file with
    the barcode conversion map if you intend to keep the same barcodes after downloading sequences from NCBI and
    demultiplexing with `harpy demultiplex ncbi`. Invalid barcodes will be generalized to 18bp of `N`, but you can
    use `--preserve-invalid`/`-p` to keep invalid barcodes unique (likely not useful for most applications).
    """
    def bx_barcode(rec):
        bx = [i for i in rec.comment.split() if i.startswith("BX:Z:")]
        #bx = re.search(r"BX:Z:[^\s]*(?=\s)", rec.comment)
        if bx:
            return bx[0].removeprefix("BX:Z:")
        else:
            return None

    def inline_barcode(rec):
        bx = re.search(r"(?:\:[ATCGN]+$|#\d+_\d+_\d+$)", rec.name)
        if bx:
            return bx[0][1:]
        else:
            return None

    def is_invalid(bx):
        return "0" in bx.split("_") or bool(re.search(r"(?:N|[ABCD]00)", bx))

    NUCLEOTIDE_FMT = False
    SPACER_NUC = "N"*5
    SPACER_QUAL = "!"*5
    ## find the barcode format
    with pysam.FastqFile(r1_fq, "r") as in_fq:
        for n,record in enumerate(in_fq, 1):
            if bx_barcode(record):
                bx_search = bx_barcode
                _bx = bx_search(record)
                if re.search(r"^[ATCGN]+$", _bx):
                    NUCLEOTIDE_FMT = True
                break
            elif inline_barcode(record):
                bx_search = inline_barcode
                _bx = bx_search(record)
                if re.search(r"^[ATCGN]+$", _bx):
                    NUCLEOTIDE_FMT = True
                break
            if n > scan:
                print_error("unknown barcode format", f"Scanned the first {scan} reads of {os.path.basename(r1_fq)} and was unable to locate barcodes in the BX:Z field nor as a TELLseq or stLFR suffix in the read ID.")

    bc_inventory = {}
    bc_iter = product(*["ATCG" for i in range(18)])
    bc_iter_inv = product(*(["N"] + ["ATCGN" for i in range(17)]))

    for i,fq in enumerate([r1_fq, r2_fq],1):
        with pysam.FastqFile(fq, "r") as in_fq, open(f"{prefix}.R{i}.fq.gz", "wb") as out_fq:
            gzip = subprocess.Popen(["gzip"], stdin = subprocess.PIPE, stdout = out_fq)
            try:
                for record in in_fq:
                    _bx = bx_search(record)
                    if not _bx:
                        inline_bc = "N"*18
                    else:
                        if NUCLEOTIDE_FMT:
                            inline_bc = _bx
                        else:
                            nuc_bx = bc_inventory.get(_bx, None)
                            if not nuc_bx:
                                if is_invalid(_bx):
                                    nuc_bx = "".join(next(bc_iter_inv)) if preserve_invalid else "N"*18
                                else:
                                    nuc_bx = "".join(next(bc_iter))
                                bc_inventory[_bx] = nuc_bx
                            inline_bc = nuc_bx
                    record.sequence = inline_bc + SPACER_NUC + record.sequence
                    record.quality  = "I"*len(inline_bc) + SPACER_QUAL + record.quality
                    gzip.stdin.write(str(record).encode("utf-8") + b"\n")
            finally:
                gzip.stdin.close()
                retcode = gzip.wait()
                if retcode != 0:
                    click.echo(f"Error: gzip exited with status {retcode}", err=True)
                    sys.exit(retcode)
    if barcode_map:
        with open(f"{prefix}.barcode.map", "w") as bc_out:
            for bx,nuc in bc_inventory.items():
                bc_out.write(f"{nuc}\t{bx}\n")

convert.add_command(fastq)
convert.add_command(ncbi)
convert.add_command(standardize_bam)
convert.add_command(standardize_fastq)
