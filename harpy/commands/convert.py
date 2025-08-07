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
from harpy.common.misc import safe_read, harpy_pulsebar
from harpy.common.convert import FQRecord, compress_fq
from harpy.common.printing import print_error
from harpy.common.validations import validate_barcodefile

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def convert():
    """
    Convert data between linked-read types
    """

module_docstring = {
    "harpy convert": [
        {
            "name": "Commands",
            "commands": ["bam", "fastq", "standardize"],
            "panel_styles": {"border_style" : "blue"}
        }
    ]
}

INVALID_10x = "N" * 16
INVALID_HAPLOTAGGING = "A00C00B00D00"
INVALID_STLFR = "0_0_0"
INVALID_TELLSEQ = "N" * 18

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/convert")
@click.option('--standardize',  is_flag = True, default = False, help = 'Add barcode validation tag `VX:i` to output')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` `1` (all) or `2` (no) output')
@click.argument('to_', metavar = 'TO', type = click.Choice(["10x","haplotagging", "stlfr", "tellseq"], case_sensitive=False), nargs = 1)
@click.argument('sam', metavar="BAM", type = click.Path(exists=True, readable=True, dir_okay=False), required = True, nargs=1)
def bam(to_,sam, standardize, quiet):
    """
    Convert between barcode formats in alignments

    This conversion changes the barcode type of the alignment file (SAM/BAM), expecting
    the barcode to be in the `BX:Z` tag of the alignment. The barcode type is automatically
    detected and the resulting barcode will be in the `BX:Z` tag. Use `--standardize` to 
    optionally standardize the output file (recommended), meaning a `VX:i` tag is added to describe
    barcode validation with `0` (invalid) and `1` (valid). Writes to `stdout`.
    """
    # Do a quick scan of the file until the first barcode
    # to assess what kind of linked read tech it is and set
    # the appropriate kind of MISSING barcode
    from_ = None
    try:
        with pysam.AlignmentFile(sam, require_index=False) as alnfile, harpy_pulsebar(quiet, "Determining barcode type", True) as progress:
            progress.add_task("[dim]Determining barcode type", total = None)
            HEADER = alnfile.header
            for record in alnfile.fetch(until_eof = True):
                try:
                    bx = record.get_tag("BX")
                    if re.search(r"^[ATCGN]+$", bx):
                        # tellseq
                        from_ = "tellseq"
                    elif re.search(r"^\d+_\d+_\d+$", bx):
                        # stlfr
                        from_ = "stlfr"
                    elif re.search(r"^A\d{2}C\d{2}B\d{2}D\d{2}$", bx):
                        from_ = "haplotagging"
                    else:
                        continue
                    break
                except KeyError:
                    continue
            if not from_:
                print_error("unrecognized barcode", f"After scanning {os.path.basename(sam)}, either no BX:Z fields were found, or no barcodes conforming to haplotagging,stlfr, or tellseq/10x were identified.")
                sys.exit(1)
    except ValueError:
        print_error("Unrecognized file type", f"[blue]{os.path.basename(sam)}[/] was unable to be processed by samtools, suggesting it is not a SAM/BAM file.")
        sys.exit(1)

    to_ = to_.lower()
    if from_ == to_:
        print_error("invalid to/from", f"The barcode formats between [green]{from_}[/] (detected from input) and [green]{to_}[/] (user-specified) must be different from each other.")
        sys.exit(1)

    # for barcodes, use sample() so the barcodes don't all start with AAAAAAAAAAAAA (or 1)
    # it's not functionally important, but it does make the barcodes *look* more distinct
    if to_ == "10x":
        bc_generator = product(*[sample("ATCG", 4) for i in range(16)])
        def is_valid(bc):
            return "N" not in bc
        def format_bc(bc):
            return "".join(bc)
    elif to_ == "tellseq":
        bc_generator = product(*[sample("ATCG", 4) for i in range(18)])
        def is_valid(bc):
            return "N" not in bc
        def format_bc(bc):
            return "".join(bc)
    elif to_ == "stlfr":
        bc_generator = product(*[sample(range(1,1538), 1537) for i in range(3)])
        def is_valid(bc):
            return "0" not in bc.split("_")
        def format_bc(bc):
            return "_".join(str(i) for i in bc)
    elif to_ == "haplotagging":
        bc_generator = product(
            ["A" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["C" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["B" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["D" + str(i).zfill(2) for i in sample(range(1,97), 96)]
        )
        def is_valid(bc):
            return "00" not in bc
        def format_bc(bc):
            return "".join(bc)

    bc_inventory = {}
    with (
        pysam.AlignmentFile(sam, require_index=False) as SAM,
        pysam.AlignmentFile(sys.stdout, "wb", header=HEADER) as OUT,
        harpy_pulsebar(quiet, "Converting", True) as progress,
    ):
        progress.add_task(f"[blue]{from_}[/] -> [magenta]{to_}[/]", total = None)
        for record in SAM.fetch(until_eof=True):
            if record.has_tag("BX"):
                bx = record.get_tag("BX")
                # the standardization is redundant but ensures being written before the BX tag
                if bx in bc_inventory:
                    if standardize:
                        record.set_tag("VX", int(is_valid(bx)), "i")
                    record.set_tag("BX", bc_inventory[bx] ,"Z")
                else:
                    try:
                        bc_inventory[bx] = format_bc(next(bc_generator))
                        if standardize:
                            record.set_tag("VX", int(is_valid(bx)), "i")
                        record.set_tag("BX", bc_inventory[bx] ,"Z")
                    except StopIteration:
                        print_error("too many barcodes", f"There are more {from_} barcodes in the input data than it is possible to generate {to_} barcodes from.")
                        sys.exit(1)
            OUT.write(record)

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/convert")
@click.option('-o','--output', type = str, metavar= "PREFIX", help='file prefix for output fastq files', required=True)
@click.option('-b','--barcodes', type = click.Path(exists=True, readable=True, dir_okay=False), help='barcodes file [from 10x only]', required=False)
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` `1` (all) or `2` (no) output')
@click.argument('from_', metavar = 'FROM', type = click.Choice(["10x", "haplotagging", "standard", "stlfr", "tellseq"], case_sensitive=False), nargs=1)
@click.argument('to_', metavar = 'TO', type = click.Choice(["10x", "haplotagging", "standard", "stlfr", "tellseq"], case_sensitive=False), nargs = 1)
@click.argument('fq1', metavar="R1_FASTQ", type = click.Path(exists=True, readable=True, dir_okay=False), required = True, nargs=1)
@click.argument('fq2', metavar="R2_FASTQ", type = click.Path(exists=True, readable=True, dir_okay=False), required=True, nargs= 1)
def fastq(from_,to_,fq1,fq2,output,barcodes, quiet):
    """
    Convert between linked-read FASTQ formats
    
    Takes the positional arguments `FROM` to indicate input data format and `TO` is the
    target data format. Both of these arguments allow the formats provided in the table below. `10x`
    input data requires a `--barcodes` file containing one nucleotide barcode per line to
    determine which barcodes are valid/invalid. In all cases, a file will be created with
    the barcode conversion map. Requires 2 threads.
    
    | from/to      | barcode format                                     | example                     |
    |:-------------|:---------------------------------------------------|:----------------------------|
    | 10x          | the first N base pairs of R1, given `--barcodes`   |                             |
    | haplotagging | a `BX:Z:ACBD` SAM tag in the sequence header       | `@SEQID BX:Z:A01C93B56D11`  |
    | standard     | a `BX:Z` SAM tag in the sequence header, any style | `@SEQID BX:Z:ATAGCAC_AGGA`  |
    | stlfr        | `#1_2_3` format appended to the sequence ID        | `@SEQID#1_2_3`              |
    | tellseq      | `:ATCG` format appended to the sequence ID         | `@SEQID:GGCAAATATCGAGAAGTC` |
    """
    if from_ == to_:
        print_error("invalid to/from", "The file formats between [green]TO[/] and [green]FROM[/] must be different from each other.")
        sys.exit(1)
    if from_ in ["haplotagging", "standard"] and to_ in ["haplotagging", "standard"]:
        print_error("redundant conversion", "The [green]haplotagging[/] and [green]standard[/] formats are functionally identical and this conversion won\'t do anything.")
        sys.exit(1)
    if from_ == "10x" and not barcodes:
        print_error("missing required file", "A [green]--barcodes[/] file must be provided if the input data is 10x.")
        sys.exit(1)
    # check that the file is fastq
    with pysam.FastxFile(fq1, persist=False) as R1:
        try:
            for i in R1:
                if i.name == "HD":
                    raise ValueError
                break
        except ValueError:
            print_error("Unrecognized file type", f"[blue]{os.path.basename(fq1)}[/] was unable to be processed as a FASTQ file by samtools, suggesting it is not a FASTQ file.")
            sys.exit(1)
    with pysam.FastxFile(fq2, persist=False) as R2:
        try:
            for i in R2:
                if i.name == "HD":
                    raise ValueError
                break
        except ValueError:
            print_error("Unrecognized file type", f"[blue]{os.path.basename(fq2)}[/] was unable to be processed as a FASTQ by samtools, suggesting it is not a FASTQ file.")
            sys.exit(1)

    # just make sure it's all lowercase
    from_ = from_.lower()
    to_ = to_.lower()
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
    elif to_ == "standard":
        if from_ == "10x":
            invalid = INVALID_10x
        elif from_ == "tellseq":
            invalid = INVALID_TELLSEQ
        elif from_ == "haplotagging":
            invalid = INVALID_HAPLOTAGGING
        elif from_ == "stlfr":
            invalid = INVALID_STLFR
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
        harpy_pulsebar(quiet, "Converting") as progress,
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
                    if (from_ in ["tellseq", "10x"] and to_ in ["tellseq", "10x"]) or to_ == "standard":
                        bc_inventory[_r1.barcode] = _r1.barcode
                    else:
                        if _r1.valid:
                            try:
                                bc_inventory[_r1.barcode] = format_bc(next(bc_generator))
                            except StopIteration:
                                print_error("too many barcodes", f"There are more {from_} barcodes in the input data than it is possible to generate {to_} barcodes from.")
                                sys.exit(1)
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
                    if (from_ in ["tellseq", "10x"] and to_ in ["tellseq", "10x"]) or to_ == "standard":
                        bc_inventory[_r2.barcode] = _r2.barcode
                    else:
                        if _r2.valid:
                            try:
                                bc_inventory[_r2.barcode] = format_bc(next(bc_generator))
                            except StopIteration:
                                print_error("too many barcodes", f"There are more {from_} barcodes in the input data than it is possible to generate {to_} barcodes from.")
                                sys.exit(1)
                        else:
                            bc_inventory[_r2.barcode] = invalid
                    bc_out.write(f"{_r2.barcode}\t{bc_inventory[_r2.barcode]}\n")
                converted_bc = bc_inventory[_r2.barcode]
                R2_out.write(str(_r2.convert(to_, converted_bc)))
    
    # bgzip compress the output, one file per thread
    with ThreadPoolExecutor(max_workers=2) as executor:
        executor.submit(compress_fq, f"{output}.R1.fq")
        executor.submit(compress_fq, f"{output}.R2.fq")

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/convert")
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` `1` (all) or `2` (no) output')
@click.argument('sam', metavar="BAM", type = click.Path(exists=True, readable=True, dir_okay=False), required = True, nargs=1)
def standardize(sam, quiet):
    """
    Move barcode to BX:Z tag in alignments

    This conversion moves the barcode from the sequence name in to the `BX:Z` tag of the alignment,
    maintaining the same barcode type (i.e. there is no linked-read format conversion). It is intended
    for tellseq and stlfr data, which encode the barcode in the read name. Also writes a `VX:i` tag
    to describe barcode validation `0` (invalid) or `1` (valid). Writes to `stdout`.
    """
    try:
        with (
            pysam.AlignmentFile(sam, require_index=False) as SAM, 
            pysam.AlignmentFile(sys.stdout, "wb", template=SAM),
            harpy_pulsebar(quiet, "Standardizing", True),
        ):
            for record in SAM.fetch(until_eof=True):
                if record.has_tag("BX") and record.has_tag("VX"):
                    print_error("BX/VX tags present", f"The BX:Z and VX:i tags are already present in {os.path.basename(sam)} and does not need to be standardized.")
                    sys.exit(1)
                if record.has_tag("BX"):
                    bx = record.get_tag("BX")
                    if "0" in bx.split("_") or re.search(r"(?:N|[ABCD]00)", bx):
                        record.set_tag("VX", 0, "i")
                    else:
                        record.set_tag("VX", 1, "i")
                    continue
                # matches either tellseq or stlfr   
                bx = re.search(r"(?:\:[ATCGN]+$|#\d+_\d+_\d+$)", record.query_name)
                if bx:
                    # the 1:0 ignores the first character, which will either be : or #
                    bx_sanitized = bx[0][1:]
                    record.query_name = record.query_name.remove_suffix(bx[0])
                    if "0" in bx_sanitized.split("_") or "N" in bx_sanitized:
                        record.set_tag("VX", 0, "i")
                    else:
                        record.set_tag("VX", 1, "i")
                    record.set_tag("BX", bx_sanitized, "Z")
    except ValueError:
        print_error("Unrecognized file type", f"[blue]{os.path.basename(sam)}[/] was unable to be processed by samtools, suggesting it is not a SAM/BAM file.")
        sys.exit(1)

@click.command(hidden = True, no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/ncbi")
@click.option('-m', '--barcode-map',  is_flag = True, default = False, help = 'Write a map of the barcode-to-nucleotide conversion')
@click.option('-p', '--prefix', required=True, type = str, help = "Output file name prefix")
@click.option('-i', '--preserve-invalid',  is_flag = True, default = False, help = 'Retain the uniqueness of invalid barcodes')
@click.option('-s', '--scan', show_default = True, default = 100, type = click.IntRange(min=1), help = 'Number of reads to scan to identify barcode location and format')
@click.argument('r1_fq', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=1)
@click.argument('r2_fq', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=1)
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
                print(f"Scanned the first {scan} reads of {os.path.basename(r1_fq)} and was unable to locate barcodes in the BX:Z field nor as a TELLseq or stLFR suffix in the read ID.")
                sys.exit(1)

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
convert.add_command(bam)
convert.add_command(ncbi)
convert.add_command(standardize)
