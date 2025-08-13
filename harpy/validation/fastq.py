
from itertools import chain
import os
import re
import pysam
from rich.markdown import Markdown
from harpy.common.printing import print_error

class FASTQ():
    '''
    A class to contain and validate FASTQ input files.
    '''
    def __init__(self, filenames):
        if any(isinstance(i, list) for i in filenames):
            self.files = list(chain.from_iterable(filenames))
        else:
            self.files = filenames
        self.bx_tag = False

        re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
        # check if any names will be clashing
        bn_r = r"[\.\_](?:[RF])?(?:[12])?(?:\_00[1-9])*?$"
        uniqs = set()
        dupes = []
        inv_pattern = r'[^a-zA-Z0-9._-]+'
        badmatch = []
        for i in self.files:
            sans_ext = os.path.basename(re_ext.sub("", str(i)))
            if sans_ext in uniqs:
                dupes.append(sans_ext)
            else:
                uniqs.add(sans_ext)
            if re.search(inv_pattern, os.path.basename(i)):
                badmatch.append(os.path.basename(i))
        if badmatch:
            print_error(
                "invalid characters",
                "Invalid characters were detected in the input FASTQ file names.",
                Markdown("Valid file names may contain only:\n- **A-Z** characters (case insensitive)\n- **.** (period)\n- **_** (underscore)\n- **-** (dash)"),
                "The offending files",
                ", ".join(badmatch)
                )
        if dupes:
            dupe_out = []
            for i in dupes:
                dupe_out.append(" ".join([j for j in self.files if i in j]))
            print_error(
                "clashing sample names",
                Markdown("Identical sample names were detected in the inputs, which will cause unexpected behavior and results.\n- files with identical names but different-cased extensions are treated as identical\n- files with the same name from different directories are also considered identical"),
                "Make sure all input files have unique names.",
                "Files with clashing names",
                dupe_out
            )
        
        self.count = len({re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in uniqs})
    
    def has_bx_tag(self, max_records: int = 50):
        """
        Parse the max_records in a list of fastq files to verify if they have BX tag (standard format). Returns as soon as the first BX tag is found.
        If a BX:Z: tag is present, updates self.bx_tag to True
        """
        for i in self.files:
            records = 0
            with pysam.FastxFile(i, persist=False) as fq:
                for record in fq:
                    records += 1
                    if "BX:Z" in record.comment:
                        self.bx_tag = True
                        return
                    if records >= max_records:
                        break
    def bc_or_bx(self, tag: str, max_records: int = 50) -> None:
        """
        Parse the list of fastq files to verify that they have BX/BC tag, and only one of those two types per file
        """
        primary = "BX:Z" if tag == "BX" else "BC:Z"
        secondary = "BC:Z" if tag == "BX" else "BX:Z"
        for fastq in self.files:
            PRIMARY = False
            SECONDARY = False
            records = 0
            with pysam.FastxFile(fastq, persist=False) as fq:
                for record in fq:
                    records += 1
                    PRIMARY = True if primary in record.comment else PRIMARY
                    SECONDARY = True if secondary in record.comment else SECONDARY
                    if PRIMARY and SECONDARY:
                        print_error(
                            "clashing barcode tags",
                            f"Both [green bold]BC:Z[/] and [green bold]BX:Z[/] tags were detected in the read headers for [blue]{os.path.basename(fastq)}[/]. Athena accepts [bold]only[/] one of [green bold]BC:Z[/] or [green bold]BX:Z[/].",
                            "Check why your data has both tags in use and remove/rename one of the tags."
                        )
                    if records >= max_records:
                        break
                # check for one or the other after parsing is done
                errtext = f" However, [green]{secondary}:Z[/] tags were detected, perhaps you meant those?" if SECONDARY else ""
                if not PRIMARY:
                    print_error(
                        "no barcodes found",
                        f"No [green bold]{primary}[/] tags were detected in read headers for [blue]{os.path.basename(fastq)}[/]. Athena requires the linked-read barcode to be present as either [green bold]BC:Z[/] or [green bold]BX:Z[/] tags.{errtext}",
                        "Check that this is linked-read data and that the barcode is demultiplexed from the sequence line into the read header as either a `BX:Z` or `BC:Z` tag."
                    )
