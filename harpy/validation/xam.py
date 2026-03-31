
from itertools import chain
import os
import pysam
import re
from harpy.common.printing import HarpyPrint

HAPLOTAGGING_RX = re.compile(r"^A\d{2}C\d{2}B\d{2}D\d{2}$")
STLFR_RX = re.compile(r"^\d+_\d+_\d+$")
TELLSEQ_RX = re.compile(r"^[ATCGN]+$")

class XAM():
    """
    A class to contain and validate BAM/SAM input files. If detect_bc is True, will scan the first `maxrec`
    records of the first [up to] 5 files to determine barcode type, stopping at the first detection of a
    recognizable barcode technology and sets the SAM.lr_type field with one of
    ["none", "haplotagging", "stlfr", "tellseq"]. The nonlinked_ok option controls whether
    the detection of "none" linked-read types is permissible, otherwise throwing an error.
    """
    def __init__(self, filenames, detect_bc:bool = False, nonlinked_ok:bool = True, check_phase:bool = False, quiet:int = 0, maxrec = 100):
        if any(isinstance(i, list) for i in filenames):
            self.files = list(chain.from_iterable(filenames))
        else:
            self.files = filenames
        self.count = 0
        self.lr_type = "none"
        self.print = HarpyPrint(quiet)
        self.max_records = maxrec
        re_ext = re.compile(r"\.(bam|sam)$", re.IGNORECASE)
        uniqs = set()
        dupes = []
        badfiles = []

        self.print.log("Validating input alignment files", newline=False)

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
            self.print.validation(False)
            self.print.error(
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
            self.print.validation(False)
            self.print.error(
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
                self.lr_type = self.which_linkedread(samfile)
                if self.lr_type != "none":
                    break
            if not nonlinked_ok and self.lr_type == "none":
                self.print.validation(False)
                self.print.error(
                    "incompatible data",
                    f"This command requires linked-read data in haplotagging, stlfr, or tellseq format, but none were found in the first {self.max_records} records of the first {len(scanned)} files.",
                    "Please double-check that these data are indeed linked-read data and the barcodes are formatted according to that technology standard.",
                    "Files Scanned",
                    "\n".join(scanned)
                )
    
        scanned: list[str] = []
        _phased: list[bool] = []
        for i,samfile in enumerate(self.files):
            if i > 5:
                break
            scanned.append(os.path.basename(samfile))
            # do lr barcode scan if enabled and not yet detected
            if detect_bc and self.lr_type == "none":
                self.lr_type = self.which_linkedread(samfile)
            # do phased scan if enabled and not yet detected
            if check_phase and not any(_phased):
                _phased.append(self.is_phased(samfile))

        if not nonlinked_ok and self.lr_type == "none":
            self.print.validation(False)
            self.print.error(
                "incompatible data",
                f"This command requires linked-read data in haplotagging, stlfr, or tellseq format, but none were found in the first {self.max_records} records of the first {len(scanned)} files.",
                "Please double-check that these data are indeed linked-read data and the barcodes are formatted according to that technology standard.",
                "Files Scanned",
                "\n".join(scanned)
            )

        if check_phase and not any(_phased):
            self.print.validation(False)
            self.print.error(
                "incompatible data",
                f"Phased alignments are required as input, but harpy was unable to find the [green]HP[/] or [green]PS[/] tags that denote phasing in the first {self.max_records} records of the first {len(scanned)} files.",
                "Please double-check that these data are indeed phased (contain [green]HP[/] or [green]PS[/] tags), otherwise you can phase these alignments using [green]harpy phase bam[/].",
                "Files Scanned",
                "\n".join(scanned)
            )
        self.print.validation(True)
    
    def is_phased(self, file_path: str) -> bool:
        ''' Scan the `file_path` to determine if the file has `PS` or `HP` tags'''
        with pysam.AlignmentFile(file_path, require_index=False) as alnfile:
            for i, record in enumerate(alnfile.fetch(until_eof = True), 1):
                if i > 100:
                    break
                if (record.has_tag("PS") or record.has_tag("HP")):
                    return True
        return False
    
    def which_linkedread(self, file_path: str) -> str:
        """
        Scans the first `self.max_records` records of a SAM/BAM file and tries to determine the barcode technology
        Returns one of: "haplotagging", "stlfr", "tellseq", or "none"
        """
        with pysam.AlignmentFile(file_path, require_index=False) as alnfile:
            for i, record in enumerate(alnfile.fetch(until_eof = True), 1):
                if i > 100:
                    break
                if not record.has_tag("BX"):
                    continue
                bx = record.get_tag("BX")
                if TELLSEQ_RX.search(bx):
                    return "tellseq"
                if STLFR_RX.search(bx):
                    return "stlfr"
                if HAPLOTAGGING_RX.search(bx):
                    return "haplotagging"
        return "none"