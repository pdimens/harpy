
from itertools import chain
import os
import pysam
import re
from harpy.common.printing import CONSOLE, print_error
from harpy.validation.barcodes import which_linkedread_sam

class XAM():
    """
    A class to contain and validate BAM/SAM input files. If detect_bc is True, will scan the first 100
    records of the first [up to] 5 files to determine barcode type, stopping at the first detection of a
    recognizable barcode technology and sets the SAM.lr_type field with one of
    ["none", "haplotagging", "stlfr", "tellseq"]. The nonlinked_ok option controls whether
    the detection of "none" linked-read types is permissible, otherwise throwing an error.
    """
    def __init__(self, filenames, detect_bc:bool = False, nonlinked_ok:bool = True, quiet:bool = False):
        if any(isinstance(i, list) for i in filenames):
            self.files = list(chain.from_iterable(filenames))
        else:
            self.files = filenames
        self.count = 0
        self.lr_type = "none"
        self.quiet: bool = quiet

        re_ext = re.compile(r"\.(bam|sam)$", re.IGNORECASE)
        uniqs = set()
        dupes = []
        badfiles = []

        if not self.quiet:
            CONSOLE.log("Validating input alignment files")

        for i in self.files:
            bn = os.path.basename(re_ext.sub("", str(i)))
            if bn in uniqs:
                dupes.append(bn)
            else:
                uniqs.add(bn)
                self.count += 1
            try:
                with pysam.AlignmentFile(i, 'r', require_index=False):
                    pass
            except (ValueError, OSError):
                badfiles.append(i)

        if badfiles:
            print_error(
                "invalid file type",
                f"[yellow]{len(badfiles)}[/] of the input alignment files did not conform to SAM/BAM format expectations.",
                "Please verify that the files listed below are properly formatted SAM or BAM files."
                "Offending Files",
                ", ".join(badfiles)
            )

        if dupes:
            dupe_out = []
            for i in dupes:
                dupe_out.append(" ".join([j for j in self.files if i in j]))
            print_error(
                "clashing sample names",
                "Identical filenames were detected, which will cause unexpected behavior and results.\n  - files with identical names but different-cased extensions are treated as identical\n  - files with the same name from different directories are also considered identical",
                "Make sure all input files have unique names.",
                "Files with Clashing Names",
                dupe_out
            )

        if detect_bc:
            scanned = []
            for i,samfile in enumerate(self.files):
                if i > 5:
                    break
                scanned.append(os.path.basename(samfile))
                self.lr_type = which_linkedread_sam(samfile)
                if self.lr_type != "none":
                    break
            if not nonlinked_ok and self.lr_type == "none":
                print_error(
                    "incompatible data",
                    "This command requires linked-read data, but harpy was unable to associate the input data as being haplotagging, stlfr, or tellseq format. Auto-detection scanned the first 100 lines of up to the first 5 files and failed to find barcodes conforming to those formatting standards.",
                    "Please double-check that these data are indeed linked-read data and the barcodes are formatted according to that technology standard.",
                    "Files Scanned",
                    "\n".join(scanned)
                )
