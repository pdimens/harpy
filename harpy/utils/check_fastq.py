import os
import re
import sys

import click
from pysam import FastxFile


class FormatChecker:
    """Base class: shared counters + SAM-spec/BX-position validation."""
    SAMSPEC = re.compile(r'[A-Z][A-Z]:[AifZHB]:')
    def __init__(self):
        self.N_READS = 0
        self.NO_BX = 0
        self.NO_VX = 0
        self.BAD_BX = 0
        self.BAD_SAM_SPEC = 0
        self.BX_NOT_LAST = 0

    def check_samspec(self, comment):
        """Validate TAG:TYPE:VALUE format and that BX:Z:, if present, is last."""
        splithead = comment.split()
        for tag in splithead:
            if not self.SAMSPEC.match(tag):
                self.BAD_SAM_SPEC += 1
        if any(tag.startswith("BX:Z:") for tag in splithead) and not splithead[-1].startswith("BX:Z:"):
            self.BX_NOT_LAST += 1

    def check_read(self, fq_record):
        """Platform-specific per-read validation. Must be overridden."""
        raise NotImplementedError

    def process(self, path):
        with FastxFile(path, persist=False) as fh:
            for entry in fh:
                self.N_READS += 1
                self.check_read(entry)
        return self

    def report(self, path):
        values = [
            os.path.basename(path),
            self.N_READS,
            self.NO_BX,
            self.NO_VX,
            self.BAD_BX,
            self.BAD_SAM_SPEC,
            self.BX_NOT_LAST,
        ]
        return "\t".join(str(v) for v in values)


class StandardChecker(FormatChecker):
    """Generic BX:Z:/VX:i: tag validation, no barcode-format check."""
    def check_read(self, fq_record):
        comment = fq_record.comment or ""
        if 'BX:Z:' not in comment:
            self.NO_BX += 1
            return
        if 'VX:i:' not in comment:
            self.NO_VX += 1
        self.check_samspec(comment)


class HaplotaggingChecker(FormatChecker):
    """Haplotagging: BX:Z: comment tag, ACBD-style barcode."""
    BARCODE_RE = re.compile(r'A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2}')
    def check_read(self, fq_record):
        comment = fq_record.comment or ""
        if 'BX:Z:' not in comment:
            self.NO_BX += 1
            return
        if 'VX:i:' not in comment:
            self.NO_VX += 1
        if not self.BARCODE_RE.search(comment):
            self.BAD_BX += 1
        self.check_samspec(comment)


class StlfrChecker(FormatChecker):
    """stLFR: barcode embedded in the read name as #n_n_n."""
    BARCODE_RE = re.compile(r'#\d+_\d+_\d+$')
    def check_read(self, fq_record):
        if not self.BARCODE_RE.search(fq_record.name):
            self.BAD_BX += 1
        self.check_samspec(fq_record.comment or "")


class TellseqChecker(FormatChecker):
    """TELL-Seq: barcode embedded in the read name as :ATCGN..."""
    BARCODE_RE = re.compile(r':[ATCGN]+$')
    def check_read(self, fq_record):
        if not self.BARCODE_RE.search(fq_record.name):
            self.BAD_BX += 1
        self.check_samspec(fq_record.comment or "")


PLATFORM_CHECKERS = {
    "standard": StandardChecker,
    "haplotagging": HaplotaggingChecker,
    "stlfr": StlfrChecker,
    "tellseq": TellseqChecker,
}


@click.command(no_args_is_help=True)
@click.argument('platform', required=True, type=click.Choice(list(PLATFORM_CHECKERS), case_sensitive=False))
@click.argument('input', required=True, type=click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden=True)
def check_fastq(platform, input):
    """
    File format validation for FASTQ file

    Parses a FASTQ file to check if any sequences don't conform to the SAM spec,
    whether BX:Z: is the last tag in the record, and the counts of: total reads,
    reads without BX:Z: tag, reads with incorrect barcode depending on the platform.
    """
    checker = PLATFORM_CHECKERS[platform.lower()]()
    checker.process(input)
    sys.stdout.write(checker.report(input) + "\n")
