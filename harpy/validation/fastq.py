
from itertools import chain
import os
import re
import pysam
from harpy.common.printing import HarpyPrint
from harpy.validation.barcodes import which_linkedread

class FASTQ():
    '''
    A class to contain and validate FASTQ input files. If detect_bc is True, will scan the first 100
    records of the first [up to] 5 forward-read files to determine barcode type, stopping at the first
    detection of a recognizable barcode technology and sets the FASTQ.lr_type field with one of
    ["none", "haplotagging", "stlfr", "tellseq"]. The nonlinked_ok option controls whether
    the detection of "none" linked-read types is permissible, otherwise throwing an error.
    '''
    def __init__(self, filenames, detect_bc:bool = False, nonlinked_ok:bool = True, quiet:int = 0):
        if any(isinstance(i, list) for i in filenames):
            self.files = list(chain.from_iterable(filenames))
        else:
            self.files = filenames
        self.bx_tag = False
        self.vx_tag = False
        # determine illumina format
        self.illumina_old = True
        self.print = HarpyPrint(quiet)
        self.lr_type = "none"
        badfiles = []

        #self.illu_old = re.compile(r'/[12]$')
        self.illu_new = re.compile(r'[12]:[YN]:\d+')
        re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
        # check if any names will be clashing
        bn_r = r"[\.\_](?:[RF])?(?:[12])?(?:\_00[1-9])*?$"
        uniqs = set()
        dupes = []
        self.print.log("Validating input FASTQ files", newline=False)

        for i in self.files:
            try:
                with pysam.FastxFile(i, persist=False) as f:
                    for j in f:
                        if not j.name or not j.quality:
                            raise ValueError
                        break
            except (ValueError, OSError):
                badfiles.append(i)

        for i in self.files:
            sans_ext = os.path.basename(re_ext.sub("", str(i)))
            if sans_ext in uniqs:
                dupes.append(sans_ext)
            else:
                uniqs.add(sans_ext)

        if badfiles:
            self.print.validation(False)
            self.print.error(
                "invalid file type",
                f"[yellow]{len(badfiles)}[/] of the input FASTQ files did not conform to format expectations.",
                "Please verify that the files listed below are properly formatted FASTQ files."
                "Offending Files",
                ", ".join(badfiles)
            )

        if dupes:
            dupe_out = []
            for i in dupes:
                dupe_out.append(" ".join([j for j in self.files if i in j]))
            self.print.validation(False)
            self.print.error(
                "clashing sample names",
                "Identical sample names were detected in the inputs, which will cause unexpected behavior and results.\n  - files with identical names but different-cased extensions are treated as identical\n  - files with the same name from different directories are also considered identical",
                "Make sure all input files have unique names.",
                "Files with Clashing Names",
                dupe_out
            )
        self.print.validation(True)
        
        self.count = len({re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in uniqs})
        self.detect_illumina_format()
        if detect_bc:
            self.print.log("Detecting linked-read barcode format", newline=False)
            scanned = []
            for i,fq in enumerate(self.files, 1):
                if i > 10:
                    break
                # skip the reverse-reads file
                if i % 2 == 0:
                    continue
                scanned.append(os.path.basename(fq))
                self.lr_type = which_linkedread(fq)
                if self.lr_type != "none":
                    break
            self.bx_tag = self.lr_type == "haplotagging"
            if not nonlinked_ok and self.lr_type == "none":
                self.print.validation(False)
                self.print.error(
                    "incompatible data",
                    "This command requires linked-read data, but harpy was unable to associate the input data as being haplotagging, stlfr, or tellseq format. Autodetection scanned the first 100 records of up to the first 5 files and failed to find barcodes conforming to those formatting standards.",
                    "Please double-check that these data are indeed linked-read data and the barcodes are formatted according to that technology standard.",
                    "Files Scanned",
                    "\n".join(scanned)
                )
            self.print.validation(True)

    def has_bx_tag(self, max_records: int = 50):
        """
        Parse the max_records in a list of fastq files to verify if they have BX tag (standard format). Returns as soon as the first BX tag is found.
        If a BX:Z: tag is present, updates self.bx_tag to True
        """
        self.print.log("Checking files for BX:Z tag", newline=False)

        for i in self.files:
            with pysam.FastxFile(i, persist=False) as fq:
                for i,record in enumerate(fq, 1):
                    if i > max_records:
                        break
                    cmt = record.comment or ""
                    if "BX:Z" in cmt:
                        self.bx_tag = True
                    if "VX:i" in cmt:
                        self.vx_tag = True
                    if self.bx_tag or self.vx_tag:
                        return
        self.print.validation(True)


    def detect_illumina_format(self):
        with pysam.FastxFile(self.files[0]) as f:
            for read in f:
                # Check old format: name ends with /1 or /2
                #if self.illu_old.search(read.name):
                #    self.illumina_old = "old"
                # Check new format: comment field starts with 1: or 2:
                if read.comment and self.illu_new.search(read.comment):
                    self.illumina_old = False
                # Only need to check the first read
                break

    def bc_or_bx(self, tag: str, max_records: int = 50) -> None:
        """
        Parse the first 50 records of a list of fastq files to verify that they have BX/BC tag, and only one of those two types per file
        """
        self.print.log("Checking files for BX:Z or BC:Z tags", newline=False)
        primary = "BX:Z" if tag == "BX" else "BC:Z"
        secondary = "BC:Z" if tag == "BX" else "BX:Z"
        for fastq in self.files:
            PRIMARY = False
            SECONDARY = False
            with pysam.FastxFile(fastq, persist=False) as fq:
                for i,record in enumerate(fq, 1):
                    if i > max_records:
                        break
                    cmt = record.comment or ""
                    PRIMARY = (primary in cmt) or PRIMARY
                    SECONDARY = (secondary in cmt) or SECONDARY
                    if PRIMARY and SECONDARY:
                        self.print.validation(False)
                        self.print.error(
                            "clashing barcode tags",
                            f"Both [green bold]BC:Z[/] and [green bold]BX:Z[/] tags were detected in the read headers for [blue]{os.path.basename(fastq)}[/]. Athena accepts [bold]only[/] one of [green bold]BC:Z[/] or [green bold]BX:Z[/].",
                            "Check why your data has both tags in use and remove/rename one of the tags."
                        )
                # check for one or the other after parsing is done
                errtext = f" However, [green]{secondary}[/] tags were detected, perhaps you meant those?" if SECONDARY else ""
                if not PRIMARY:
                    self.print.validation(False)
                    self.print.error(
                        "no barcodes found",
                        f"No [green bold]{primary}[/] tags were detected in read headers for [blue]{os.path.basename(fastq)}[/]. Athena requires the linked-read barcode to be present as either [green bold]BC:Z[/] or [green bold]BX:Z[/] tags.{errtext}",
                        "Check that this is linked-read data and that the barcode is demultiplexed from the sequence line into the read header as either a `BX:Z` or `BC:Z` tag."
                    )
        self.print.validation(True)
