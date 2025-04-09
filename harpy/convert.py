#! /usr/bin/env python3

from concurrent.futures import ThreadPoolExecutor
import gzip
from itertools import product, zip_longest
import os
from random import sample
import sys
import rich_click as click
import pysam
from ._misc import safe_read, harpy_pulsebar
from ._cli_types_generic import convert_to_int
from ._printing import print_error
from ._validations import validate_barcodefile

INVALID_10x = "N" * 16
INVALID_HAPLOTAGGING = "A00C00B00D00"
INVALID_STLFR = "0_0_0"
INVALID_TELLSEQ = "N" * 18

#TODO figure out how to strip and then add the 1:N:blahblah
class FQRecord():
    def __init__(self, pysamfq, FORWARD: bool, bc: str, length: int):
        self.forward = FORWARD
        fr = "1:N:" if self.forward else "2:N:"
        comments = pysamfq.comment.strip().split()
        designation = [i for i in comments if i.startswith(fr)]
        self.illumina_new = designation[0] if designation else f"{fr}0:CAGATC"
        self.illumina_old = "/1" if self.forward else "/2"
        self.id = pysamfq.name.rstrip(self.illumina_old)
        self.comment = "\t".join(i for i in comments if not i.startswith("1:N:"))
        self.seq = pysamfq.sequence
        self.qual = pysamfq.quality
        self.valid = True
        # account for /1 and /2 formatting at the end of the read (if present)
        if bc == "stlfr":
            # identify the trailing #1_2_3 barcode and remove it from the ID
            if self.id.endswith("#"):
                # is invalid
                self.id = self.id.rstrip("#")
                self.barcode = INVALID_STLFR
                self.valid = False
            else:
                _id = pysamfq.name.split("#")
                self.barcode = _id.pop(-1)
                self.id = "#".join(_id)
                self.valid = "_0" not in self.barcode
        elif bc == "10x":
            # identify the first N bases and remove it from the seq and qual of R1
            if self.forward:
                self.seq = pysamfq.sequence[length:]
                self.qual = pysamfq.quality[length:]
                self.barcode = pysamfq.sequence[:length]
        elif bc == "tellseq":
            # identify the trailing :ATCG barcode and remove it from the ID
            if self.id.endswith(":"):
                # is invalid
                self.id = self.id.rstrip(":")
                self.barcode = INVALID_TELLSEQ
                self.valid = False
            else:
                _id = pysamfq.name.split(":")
                self.barcode = _id.pop(-1)
                self.id = ":".join(_id)
                self.valid = "N" not in self.barcode
        elif bc in ["haplotagging", "standard"]:
            # identify the BX:Z SAM tag and remove it from the comment
            bc = [i for i in self.comment.split() if i.startswith("BX:Z")]
            if len(bc) < 1:
                # is invalid
                self.barcode = INVALID_HAPLOTAGGING
            else:
                self.barcode = bc.pop().removeprefix("BX:Z:")
            self.comment = "\t".join(i for i in self.comment.split() if not i.startswith("BX:Z"))
            self.valid = "00" not in self.barcode
        else:
            # the barcode was provided outright, so just use it
            # this is used in cases where the R2/R1 doesn't have retrievable barcode info (10X R2)
            self.barcode = bc
            # clear out existing BX or # tags
            self.comment = "\t".join(i for i in self.comment.split() if not i.startswith("BX:Z"))
            if "#" in self.id:
                _id = pysamfq.name.split("#")
                bc = _id.pop(-1)
                self.id = "#".join(_id)
            self.id = self.id.rstrip(":")
    
    def __str__(self):
        """Default string method returns a formatted FASTQ record."""
        return f"@{self.id}\t{self.comment}\n{self.seq}\n+\n{self.qual}\n"

    def convert(self, _type: str, BC: str):
        if _type == "10x":
            if self.forward:
                self.seq = BC + self.seq
                self.qual = "F"*len(BC) + self.qual
            self.id += f" {self.illumina_new}"
        elif _type == "tellseq":
            self.id += f":{BC} {self.illumina_new}"
        elif _type == "stlfr":
            self.id += f"#{BC} {self.illumina_new}"
        elif _type in ["haplotagging", "standard"]:
            self.id += self.illumina_old
            self.comment += f"\tBX:Z:{BC}"
        return self

def compress_fq(fq: str):
    """use pysam bgzip to compress fastq and delete the original"""
    pysam.tabix_compress(fq, f"{fq}.gz", force=True)
    os.remove(fq)

@click.command(epilog = "Documentation: https://pdimens.github.io/harpy/convert")
@click.option('-o','--output', type = str, metavar= "PREFIX", help='file prefix for output fastq files', required=True)
@click.option('-b','--barcodes', type = click.Path(exists=True, readable=True, dir_okay=False), help='barcodes file [from 10x only]', required=False)
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "2"]), callback = convert_to_int, help = '`0` (all) or `2` (no) output')
@click.argument('from_', metavar = 'FROM', type = click.Choice(["10x", "haplotagging", "standard", "stlfr", "tellseq"], case_sensitive=False), nargs=1)
@click.argument('to_', metavar = 'TO', type = click.Choice(["10x", "haplotagging", "standard", "stlfr", "tellseq"], case_sensitive=False), nargs = 1)
@click.argument('fq1', metavar="R1_FASTQ", type = click.Path(exists=True, readable=True, dir_okay=False), required = True, nargs=1)
@click.argument('fq2', metavar="R2_FASTQ", type = click.Path(exists=True, readable=True, dir_okay=False), required=True, nargs= 1)
def convert(from_,to_,fq1,fq2,output,barcodes, quiet):
    """
    Convert between linked-read FASTQ formats
    
    Takes the positional arguments `FROM` to indicate input data format and `TO` is the
    target data format. Both of these arguments allow the formats provided in the table below. `10x`
    input data requires a `--barcodes` file containing one nucleotide barcode per line to
    determine which barcodes are valid/invalid. In all cases, a file will be created with
    the barcode conversion map. Requires 2 threads.
    
    | from/to | barcode format | example |
    |:------|:-------|:--------|
    |10x    |the first N base pairs of R1, given `--barcodes` | |
    |haplotagging | a `BX:Z:ACBD` SAM tag in the sequence header | `@SEQID BX:Z:A01C93B56D11` |
    | standard | a `BX:Z` SAM tag in the sequence header, any style | `@SEQID BX:Z:ATAGCAC_AGGA` |
    | stlfr | `#1_2_3` format appended to the sequence ID | `@SEQID#1_2_3` |
    | tellseq | `:ATCG` format appended to the sequence ID | `@SEQID:GGCAAATATCGAGAAGTC` |
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
                _r2 = FQRecord(r2, False, _r1.barcode if from_ == "10x" else from_, bc_len)
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
