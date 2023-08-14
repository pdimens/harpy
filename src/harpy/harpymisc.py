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

def createregions(infile, window, base):
    bn = os.path.basename(infile)
    os.makedirs("Genome", exist_ok = True)
    if not os.path.exists(f"Genome/{bn}"):
        shell(f"ln -sr {infile} Genome/{bn}")
    if not os.path.exists(f"Genome/{bn}.fai"):
        print(f"Genome/{bn}.fai not found, indexing {bn} with samtools faidx", file = sys.stderr)
        subprocess.run(["samtools","faidx", "--fai-idx", f"Genome/{bn}.fai", infile, "2>", "/dev/null"])
    with open(f"Genome/{bn}.fai") as fai:
        bedregion = []
        while True:
            # Get next line from file
            line = fai.readline()
            # if line is empty, end of file is reached
            if not line:
                break
            # split the line by tabs
            lsplit = line.split()
            contig = lsplit[0]
            c_len = int(lsplit[1])
            start = base
            end = window
            starts = [base]
            ends = [window]
            while end < c_len:
                end = end + window if (end + window) < c_len else c_len
                ends.append(end)
                start += window
                starts.append(start)
            for (startpos, endpos) in zip (starts,ends):
                bedregion.append(f"{contig}:{startpos}-{endpos}")
        return bedregion

_regions   = createregions(genomefile, chunksize)

def check_impute_params(parameters):
    with open(parameters, "r") as fp:
        header = fp.readline().rstrip().lower()
        headersplt = header.split()
        correct_header = sorted(["model", "usebx", "k", "s", "ngen"])
        row = 0
        badrows = []
        badlens = []
        if sorted(headersplt) != correct_header:
            culprits = [i for i in headersplt if i not in correct_header]
            print(f"\n\033[1;33mERROR:\033[00m Parameter file \033[01m{parameters}\033[00m has incorrect column names. Valid names are:\n\tmodel usebx k s ngen\n", file = sys.stderr)
            print("Column names causing this error:\n\t" + " ".join(culprits), file = sys.stderr)
            print(f"\n\033[1;34mSOLUTION:\033[00m Fix the headers in \033[01m{parameters}\033[00m or use \033[01mharpy extra -s stitch.params\033[00m to generate a valid parameter file and modify it with appropriate values.")
            sys.exit(1)
        # instantiate dict with colnames
        data = dict()
        for i in headersplt:
            data[i] = []
        while True:
            # Get next line from file
            line = fp.readline()
            row += 1
            # if line is empty, end of file is reached
            if not line:
                break
            # split the line by whitespace
            rowvals = line.rstrip().split()
            rowlen = len(rowvals)
            if rowlen != 5:
                badrows.append(row)
                badlens.append(rowlen)
            elif len(badrows) > 0:
                continue
            else:
                for j,k in zip(headersplt, rowvals):
                    data[j].append(k)

        if len(badrows) > 0:
            print(f"\n\033[1;33mERROR:\033[00m Parameter file \033[01m{parameters}\033[00m is formatted incorrectly. Not all rows have the expected 5 columns.", file = sys.stderr)
            print(f"\n\033[1;34mSOLUTION:\033[00m See the problematic rows below. Check that you are using a whitespace (space or tab) delimeter in \033[01m{parameters}\033[00m or use \033[01mharpy extra -s stitch.params\033[00m to generate a valid parameter file and modify it with appropriate values.")
            print("\033[01mrow\tcolumns\033[00m")
            for i in zip(badrows, badlens):
                print(f"{i[0]}\t{i[1]}")
            sys.exit(1)
        
        # Validate each column
        culprits = {
            "model" : [],
            "usebx" : [],
            "k"     : [],
            "s"     : [],
            "ngen"  : []
        }
        colerr = 0
        for i,j in enumerate(data["model"]):
            if j not in ["pseudoHaploid", "diploid","diploid-inbred"]:
                culprits["model"].append(str(i + 1))
                colerr += 1
        if culprits["model"]:
            print("Invalid values for column \033[01mmodel\033[00m.", file = sys.stderr)
            print("Expected values: diploid, diploid-inbred, pseudoHaploid (case-sensitive)", file = sys.stderr)
            print("Rows causing error: " + " ".join(culprits["model"]), file = sys.stderr)
            print("", file = sys.stderr)

        for i,j in enumerate(data["usebx"]):
            if j not in [True, "TRUE", "true", False, "FALSE", "false", "Y","y", "YES", "Yes", "yes", "N", "n", "NO", "No", "no"]:
                culprits["usebx"].append(str(i + 1))
                colerr += 1
        if culprits["usebx"]:
            print("Invalid values for column \033[01musebx\033[00m.", file = sys.stderr)
            print("Expected values: True, False (not case sensitive)", file = sys.stderr)
            print("Rows causing error: " + " ".join(culprits["usebx"]), file = sys.stderr)
            print("", file = sys.stderr)

        for i,j in enumerate(data["k"]):
            if not j.isdigit():
                culprits["k"].append(str(i + 1))
                colerr += 1

        if culprits["k"]:
            print("Invalid values for column \033[01mk\033[00m.", file = sys.stderr)
            print("Expected values: Integers", file = sys.stderr)
            print("Rows causing error: " + " ".join(culprits["k"]), file = sys.stderr)
            print("", file = sys.stderr)

        for i,j in enumerate(data["s"]):
            if not j.isdigit():
                culprits["s"].append(str(i + 1))
                colerr += 1
        if culprits["s"]:
            print("Invalid values for column \033[01ms\033[00m.", file = sys.stderr)
            print("Expected values: Integers", file = sys.stderr)
            print("Rows causing error: " + " ".join(culprits["s"]), file = sys.stderr)
            print("", file = sys.stderr)

        for i,j in enumerate(data["ngen"]):
            if not j.isdigit():
                culprits["ngen"].append(str(i + 1))
                colerr += 1
        if culprits["ngen"]:
            print("Invalid values for column \033[01mngen\033[00m.", file = sys.stderr)
            print("Expected values: Integers", file = sys.stderr)
            print("Rows causing error: " + " ".join(culprits["ngen"]), file = sys.stderr)
            print("", file = sys.stderr)

        if colerr > 0:
            print(f"\033[1;33mERROR:\033[00m Parameter file \033[01m{parameters}\033[00m is formatted incorrectly. Not all columns have valid values.", file = sys.stderr)
            print(f"\n\033[1;34mSOLUTION:\033[00m See above for an explanation of what values each column expects and which rows are causing problems.")
            print("Find more details in the documentation: https://pdimens.github.io/harpy/modules/impute/#parameter-file", file = sys.stderr)
            exit(1)
        
# standardizes fastq filenames to not include periods in the samplename and
# ends in .F|R.fq[.gz]
# dont use?
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
        except:
            pass
    return samplenames