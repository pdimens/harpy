
import os
import re

import pysam

from harpy.common.printing import HarpyPrint

alphanum = re.compile(r'^[a-zA-Z0-9_\-\|]+$')
nuc = re.compile(r'[ACGTURYKMSWBDHVN\-]+$', flags = re.IGNORECASE)

class FASTA():
    '''
    A class to contain and validate a FASTA input file.
    '''
    def __init__(self, fasta, quiet:int = 0):
        self.file = fasta
        self.file_base = os.path.basename(self.file)
        self.print = HarpyPrint(quiet)
        self.print.log("Validating input FASTA file", newline=False)
        try:
            with pysam.FastxFile(self.file, persist=False) as fa:
                pass
        except Exception as e:
            self.print.validation(False)
            self.print.error("bad FASTA file", e)
        self.print.validation(True)


    def match_contigs(self, contigs: str):
        """Checks whether a list of contigs are present in a fasta file"""
        self.print.log("Contigs present in input FASTA", newline=False)
        valid_contigs = []
        with pysam.FastxFile(self.file, persist=False) as fa:
            for record in fa:
                valid_contigs.append(record.name)
        bad_names = []
        for i in contigs:
            if i not in valid_contigs:
                bad_names.append(i)
        if bad_names:
            self.print.validation(False)
            self.print.error(
                "contigs absent",
                f"Some of the provided contigs were not found in [blue]{self.file_base}[/]. This will definitely cause plotting errors in the workflow.",
                "Check that your contig names are correct, including uppercase and lowercase.",
                f"Contigs absent in {self.file_base}",
                ",".join([i for i in bad_names])
            )
        self.print.validation(True)

    def contigs(self) -> dict:
        """Read the FASTA file to retrieve contig names and lengths"""
        contigs = {}
        with pysam.FastxFile(self.file, persist=False) as fopen:
            for record in fopen:
                contigs[record.name] = len(record.sequence)
        return contigs

    def validate_region(self, regioninput) -> None:
        """validates the --regions input of harpy snp to infer whether it's an integer, region, or file"""
        try:
            # is an int
            int(regioninput)
            return
        except ValueError:
            pass

        self.print.log("Input regions present in FASTA", newline=False)

        # is a file specifying regions
        contigs = self.contigs()
        if os.path.isfile(regioninput):
            with open(regioninput, "r", encoding="utf-8") as fin:
                for idx, line in enumerate(fin, 1):
                    row = line.split()
                    if len(row) != 3:
                        self.print.validation(False)
                        self.print.error(
                            "invalid format",
                            f"The input file is formatted incorrectly at line {idx}. This is the first row triggering this error, but it may not be the only one.",
                            f"Rows in [blue]{regioninput}[/] need to be [bold]space[/] or [bold]tab[/] delimited with the format [yellow bold]contig start end[/] where [yellow bold]start[/] and [yellow bold]end[/] are integers.",
                            "Rows triggering this error",
                            line
                        )
                    else:
                        try:
                            start = int(row[1])
                            end = int(row[2])
                        except ValueError:
                            self.print.error(
                                "invalid format",
                                f"The input file is formatted incorrectly at line {idx}. This is the first row triggering this error, but it may not be the only one.",
                                f"Rows in [blue]{regioninput}[/] need to be [bold]space[/] or [bold]tab[/] delimited with the format [yellow bold]contig start end[/] where [yellow bold]start[/] and [yellow bold]end[/] are integers.",
                                "Rows triggering this error",
                                line
                            )
                    if row[0] not in contigs:
                        self.print.error(
                            "missing contig",
                            f"The contig listed at row {idx} ([bold yellow]{row[0]}[/]) is not present in ([blue]{self.file_base}[/]). This is the first row triggering this error, but it may not be the only one.",
                            f"Check that all the contigs listed in [blue]{os.path.basename(regioninput)}[/] are also present in [blue]{self.file_base}[/]",
                            "Row triggering this error",
                            line
                        )
                    if start > end:
                        self.print.error(
                            "invalid interval",
                            f"The interval start position is greater than the interval end position at row {idx}. This is the first row triggering this error, but it may not be the only one.",
                            f"Check that all rows in [blue]{os.path.basename(regioninput)}[/] have a [bold yellow]start[/] position that is less than the [bold yellow]end[/] position.",
                            "Row triggering this error",
                            line
                        )
                    if start > contigs[row[0]] or end > contigs[row[0]]:
                        self.print.error(
                            "invalid interval",
                            f"The interval start or end position is out of bounds at row {idx}. This is the first row triggering this error, but it may not be the only one.",
                            f"Check that the intervals present in [blue]{os.path.basename(regioninput)}[/] are within the bounds of the lengths of their respective contigs. This specific error is triggered for [bold yellow]{row[0]}[/], which has a total length of [bold]{contigs[row[0]]}[/].",
                            "Row triggering this error",
                            line
                        )
            self.print.validation(True)
            return

        contig,positions = regioninput.split(":")
        startpos,endpos = [int(i) for i in positions.split("-")]
        if contig not in contigs:
            self.print.validation(False)
            self.print.error("contig not found", f"The contig ([bold yellow]{contig}[/]) of the input region [yellow bold]{regioninput}[/] was not found in [blue]{self.file}[/].")
        if startpos > contigs[contig]:
            self.print.validation(False)
            self.print.error("region out of bounds", f"The start position ([yellow bold]{startpos}[/]) exceeds the total length of contig [yellow bold]{contig}[/] ({contigs[contig]})")
        if endpos > contigs[contig]:
            self.print.validation(False)
            self.print.error("region out of bounds", f"The end position ([yellow bold]{endpos}[/]) exceeds the total length of contig [yellow bold]{contig}[/] ({contigs[contig]})")

        self.print.validation(True)
        return
