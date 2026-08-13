import re
import sys

import click
from pysam import FastxFile

class haplotagging:
    def __init__(self):
        self.bc = re.compile(r'BX:Z:(A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2})')
        self.invalids: dict[str, int] = {"A" : 0, "C" : 0, "B" : 0, "D" : 0}
        self.bcInvalid = re.compile('([AaBbCcDd])00')
    def process_invalid(self, barcode) -> bool:
        for i in self.bcInvalid.findall(barcode):
            self.invalids[i] += 1
        return '00' in barcode
    def parse(self, read):
        '''search for a haplotagging barode in the read comment'''
        return self.bc.search(read.comment)

class stlfr:
    def __init__(self):
        self.bc = re.compile('#([0-9]+_[0-9]+_[0-9]+)')
        self.invalids: dict[str, int] = {"1" : 0, "2" : 0, "3" : 0}
    def process_invalid(self, barcode):
            split_bc = barcode.split("_")
            for i,j in enumerate(split_bc, 1):
                if j == "0":
                    self.invalids[str(i)] += 1
            return '0' in split_bc
    def parse(self, read):
        ''' search for a stlfr barode in the read name'''
        return self.bc.search(read.name)

class tellseq:
    def __init__(self):
        self.bc = re.compile(':([ATCGN]+)')
        self.invalids: dict[str, int] = dict()
        for i in range(1,19):
            self.invalids[str(i)] = 0
        self.bcInvalid = re.compile('N', flags=re.IGNORECASE)
    def process_invalid(self, barcode):
        for i in self.bcInvalid.finditer(barcode):
            self.invalids[str(i.start()+1)] += 1
        return 'N' in barcode
    def parse(self, read):
        ''' search for a tellseq barode in the read name'''
        return self.bc.search(read.name)

@click.command(no_args_is_help = True)
@click.argument('platform', required = True, type=click.Choice(['haplotagging','stlfr','tellseq'], case_sensitive=False))
@click.argument('input', required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def bx_stats_fq(platform, input):
    """
    Parses a FASTQ file to count: total sequences, total number of linked-read barcodes,
    number of valid barcodes, number of invalid BX tags, and a count of positional
    barcode invalidations (e.g. A00, _0_, N)
    """
    if platform == "haplotagging":
        parser = haplotagging()
    elif platform == "stlfr":
        parser = stlfr()
    else:
        parser = tellseq()

    N_READS = 0
    N_BX = 0
    N_VALID = 0

    with FastxFile(input, persist = False) as fh:
        for entry in fh:
            N_READS += 1
            query = parser.parse(entry)
            if query:
                N_BX += 1
                barcode = query.group(1)
                is_invalid = parser.process_invalid(barcode)
                if not is_invalid:
                    N_VALID += 1

    sys.stdout.write(f"Reads\t{N_READS}\n")
    sys.stdout.write(f"Barcodes\t{N_BX}\n")
    sys.stdout.write(f"Barcodes_Valid\t{N_VALID}\n")
    sys.stdout.write(f"Barcodes_Invalid\t{N_BX - N_VALID}\n")
    for bx,count in parser.invalids.items():
        sys.stdout.write(f"{bx}\t{count}\n")
