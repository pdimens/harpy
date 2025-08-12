"""custom classes and functions for convert"""

import os
import pysam
from .printing import print_error

INVALID_STLFR = "0_0_0"
INVALID_TELLSEQ = "N" * 18
INVALID_HAPLOTAGGING = "A00C00B00D00"

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
    try:
        pysam.tabix_compress(fq, f"{fq}.gz", force=True)
        os.remove(fq)
    except Exception as e:
        print_error("compression error", f"Failed to compress {fq}: {str(e)}")
