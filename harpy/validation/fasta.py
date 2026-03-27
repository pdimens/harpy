
import os
import pysam
from harpy.common.file_ops import safe_read
from harpy.common.printing import HarpyPrint

class FASTA():
    '''
    A class to contain and validate a FASTA input file.
    '''
    def __init__(self, fasta, quiet:int = 0):
        self.file = fasta
        self.print = HarpyPrint(quiet)
        self.print.log("Validating input FASTA file format", newline = False)
        # validate fasta file contents
        line_num = 0
        seq_id = 0
        seq = 0
        last_header = False
        with safe_read(self.file) as fasta:
            for line in fasta:
                line_num += 1
                if line.startswith(">"):
                    seq_id += 1
                    if last_header:
                        self.print.validation(False)
                        self.print.error(
                            "consecutive contig names",
                            f"All contig names must be followed by at least one line of nucleotide sequences, but two consecutive lines of contig names were detected. This issue was identified at line [bold]{line_num}[/] in [blue]{self.file}[/], but there may be others further in the file.",
                            "See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/"
                        )
                    else:
                        last_header = True
                    if len(line.rstrip()) == 1:
                        self.print.validation(False)
                        self.print.error(
                            "unnamed contigs",
                            f"All contigs must have an alphanumeric name, but a contig was detected without a name. This issue was identified at line [bold]{line_num}[/] in [blue]{self.file}[/], but there may be others further in the file.",
                            "See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/"
                        )
                    if line.startswith("> "):
                        self.print.validation(False)
                        self.print.error(
                            "invalid contig names",
                            f"All contig names must be named [green bold]>contig_name[/], without a space, but a contig was detected with a space between the [green bold]>[/] and contig_name. This issue was identified at line [bold]{line_num}[/] in [blue]{self.file}[/], but there may be others further in the file.",
                            "See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/"
                        )
                elif line == "\n":
                    self.print.validation(False)
                    self.print.error(
                        "empty lines",
                        f"Empty lines are not permitted in FASTA files, but one was detected at line [bold]{line_num}[/] in [blue]{self.file}[/]. The scan ended at this error, but there may be others further in the file.",
                        "See the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/"
                    )
                else:
                    seq += 1
                    last_header = False
        solutiontext = "FASTA files must have at least one contig name followed by sequence data on the next line. Example:\n"
        solutiontext += "[green]  >contig_name\n  ATACAGGAGATTAGGCA[/]\n"
        # make sure there is at least one of each
        if seq_id == 0:
            self.print.validation(False)
            self.print.error(
                "contig names absent",
                f"No contig names detected in [blue]{self.file}[/].",
                f"{solutiontext}\nSee the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/"
            )
        if seq == 0:
            self.print.validation(False)
            self.print.error(
                "sequences absent",
                f"No sequences detected in [blue]{self.file}[/].",
                f"{solutiontext}\nSee the FASTA file spec and try again after making the appropriate changes: https://www.ncbi.nlm.nih.gov/genbank/fastaformat/"
            )
        self.print.validation(True)


    def match_contigs(self, contigs: str):
        """Checks whether a list of contigs are present in a fasta file"""
        self.print.log("Validating input contigs against those in the input FASTA", newline=False)
        valid_contigs = []
        with pysam.FastxFile(self.file, persist=False) as fa:
            for record in fa:
                valid_contigs.append(record.name)
        bad_names = []
        for i in contigs:
            if i not in valid_contigs:
                bad_names.append(i)
        if bad_names:
            shortname = os.path.basename(self.file)
            self.print.validation(False)
            self.print.error(
                "contigs absent",
                f"Some of the provided contigs were not found in [blue]{shortname}[/]. This will definitely cause plotting errors in the workflow.",
                "Check that your contig names are correct, including uppercase and lowercase.",
                f"Contigs absent in {shortname}",
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

        self.print.log("Validating input regions to those in the input FASTA", newline=False)

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
                            f"The contig listed at row {idx} ([bold yellow]{row[0]}[/]) is not present in ([blue]{os.path.basename(self.file)}[/]). This is the first row triggering this error, but it may not be the only one.",
                            f"Check that all the contigs listed in [blue]{os.path.basename(regioninput)}[/] are also present in [blue]{os.path.basename(self.file)}[/]",
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