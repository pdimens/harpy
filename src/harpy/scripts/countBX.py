import pysam
import re
import sys

n_reads = 0
n_bx = 0
n_valid = 0
# haplotag = re.compile("([A-Z]\d{2,}){3,}")
haplotag = re.compile('A[0-9]{2}C[0-9]{2}B[0-9]{2}D[0-9]{2}')
# invalid = re.compile('[A-Z]00')
invalid = re.compile('[AaBbCcDd]00')
# inv_dict = dict()
inv_dict = {
    "A" : 0,
    "B" : 0,
    "C" : 0,
    "D" : 0
}
with pysam.FastxFile(snakemake.input[0]) as fh:
    for entry in fh:
        n_reads += 1
        comments = entry.comment.split()
        # looking for a comment that starts as whitespace + BX:Z:
        bxtag_idx = [i for i,j in enumerate(comments) if j.startswith("BX:Z:")]
        #if 'BX:Z:' in entry.comment:
        if bxtag_idx:
            n_bx += 1
            beadtag_full = comments[bxtag_idx[0]]
            beadtag = beadtag_full[5:]
            if bool(haplotag.match(beadtag)):
                inv = re.findall(invalid, beadtag)
                if inv:
                    for i in inv:
                    #    if i[0] in inv_dict:
                        inv_dict[i[0]] += 1
                    #   else:
                    #   inv_dict[i[0]] = 1
                    continue
                n_valid += 1

with open(snakemake.output[0], "w") as fout:
    print(f"totalReads\t{n_reads}", file = fout)
    print(f"bxTagCount\t{n_bx}", file = fout)
    print(f"bxValid\t{n_valid}", file = fout)
    print(f"bxInvalid\t{n_bx - n_valid}", file = fout)
    print("A00\t",str(inv_dict["A"]), file = fout)
    print("C00\t",str(inv_dict["C"]), file = fout)
    print("B00\t",str(inv_dict["B"]), file = fout)
    print("D00\t",str(inv_dict["D"]), file = fout)