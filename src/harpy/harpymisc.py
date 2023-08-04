import sys
import os

## recurring checks and such ##
def vcfcheck(vcf):
    vfile = vcf.lower()
    if vfile.endswith(".vcf") or vfile.endswith(".bcf") or vfile.endswith(".vcf.gz"):
        pass
    else:
        print(f"\033[1;33mERROR:\033[00m Supplied variant call file ({vcf}) must end in one of [.vcf | .vcf.gz | .bcf]", file = sys.stderr)
        exit(1)

def getnames(directory, ext):
    samplenames = set([i.split(ext)[0] for i in os.listdir(directory) if i.endswith(ext)])
    if len(samplenames) < 1:
        print(f"\033[1;33mERROR:\033[00m No sample files ending with {ext} found in {directory}.", file = sys.stderr)
        sys.exit(1)
    return samplenames

# Nicer version for init
def getnames_err(directory, ext):
    samplenames = set([i.split(ext)[0] for i in os.listdir(directory) if i.endswith(ext)])
    if len(samplenames) < 1:
        sys.tracebacklimit = 0
        raise Exception(f"\033[1;33mERROR:\033[00m No sample files ending with {ext} found in {directory}.")
    return samplenames
