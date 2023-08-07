import sys
import os
import re
import gzip
from pathlib import Path

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

# standardizes fastq filenames to not include periods in the samplename and
# ends in .F|R.fq[.gz]
def sanitize_fastq(full_fqlist, linkdir):
    samplenames = set()
    os.makedirs(linkdir, exist_ok = True)
    # regex to find forward or reverse spec in read header
    bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
    re_FR = re.compile(r"\s[12]\:[YN]")
    for seqfile in full_fqlist:
        if seqfile.endswith(".gz"):
            gz_ext = ".gz"
            f = gzip.open(seqfile, "r")
        else:
            gz_ext = ""
            f = open(seqfile, "r")
        header = f"{f.readline()}"
        f.close()
        # search for 1:Y or 2:Y in first read header
        if "1" in re_FR.search(header).group(0):
            FR = "F"
        else:
            FR = "R"
        bn = re.sub(bn_r, "", os.path.basename(seqfile), flags = re.IGNORECASE)
        bn = bn.replace(".", "_")
        samplenames.add(bn)
        target = Path(seqfile).absolute()
        linkedfile = f"{linkdir}/{bn}.{FR}.fq{gz_ext}"
        try:
            _ = Path(linkedfile).symlink_to(target)
            #subprocess.run(["ln", "-sr", seqfile, f"{linkdir}/{bn}.{FR}.fq{gz_ext}"])
        except:
            pass
    return samplenames