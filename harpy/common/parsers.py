import re

class LRparser:
    """Base class for FASTQ/SAM parsing"""

class Haplotagging(LRparser):
    '''Haplotagging data parser'''
    def __init__(self):
        self.patternFQ = re.compile(r'BX:Z:(A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2})')
        self.patternInvalid = re.compile('([AaBbCcDd])00')
        self.invalidDict: dict[str, int] = {"A" : 0, "C" : 0, "B" : 0, "D" : 0}
    def process_invalid(self, barcode) -> bool:
        for i in self.patternInvalid.findall(barcode):
            self.invalidDict[i] += 1
        return '00' in barcode
    def getBX_fq(self, read):
        '''search for a haplotagging barode in the read comment'''
        return self.patternFQ.search(read.comment)

class Stlfr(LRparser):
    '''stLFR data parser'''
    def __init__(self):
        self.patternFQ = re.compile('#([0-9]+_[0-9]+_[0-9]+)')
        self.invalidDict: dict[str, int] = {"1" : 0, "2" : 0, "3" : 0}
    def process_invalid(self, barcode):
            split_bc = barcode.split("_")
            for i,j in enumerate(split_bc, 1):
                if j == "0":
                    self.invalidDict[str(i)] += 1
            return '0' in split_bc
    def getBX_fq(self, read):
        ''' search for a stlfr barode in the read name'''
        return self.patternFQ.search(read.name)

class Tellseq(LRparser):
    '''TELLseq data parser'''
    def __init__(self):
        self.patternFQ = re.compile(':([ATCGN]+)')
        self.invalidDict: dict[str, int] = dict()
        for i in range(1,19):
            self.invalidDict[str(i)] = 0
        self.patternInvalid = re.compile('N', flags=re.IGNORECASE)
    def process_invalid(self, barcode):
        for i in self.patternInvalid.finditer(barcode):
            self.invalidDict[str(i.start()+1)] += 1
        return 'N' in barcode
    def getBX_fq(self, read):
        ''' search for a tellseq barode in the read name'''
        return self.patternFQ.search(read.name)