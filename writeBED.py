import subprocess

max_distance = 50000
outfile = str(input().strip())
outfile = outfile.replace(".bam", ".BX.bed")
with open(outfile, "w") as OUT:
    last_bx = ""
    last_read = ""
    last_chr = ""
    cache = []
    sort_order = {}
    with subprocess.Popen(["samtools", "view", outfile, "-F", "1024"], stdout=subprocess.PIPE) as pipe:
        for line in pipe.stdout:
            line = line.decode().strip()
            temp = line.split("\t")
            bx = re.search(r"BX:Z:(\S+)", line).group(1)
            if temp[5] == "*" or bx.find("A00" * 3) != -1 or bx.find("A00" * 4) != -1:
                continue
            if bx != last_bx or temp[2] != last_chr:
                mols = []
                last_end = 0
                if cache:
                    str = "\t".join(cache[0])
                    for ci in range(1, len(cache)):
                        str = "\t".join(cache[ci])
                        last_cigar = re.findall(r"(\d+[MIDNSHPX=])", cache[ci-1][5])
                        last_a_len = 0
                        for v in last_cigar:
                            if re.match(r"(\d+)[MDNX=]", v):
                                last_a_len += int(re.match(r"(\d+)[MDNX=]", v).group(1))
                        if cache[ci][3] - (cache[ci-1][3] + last_a_len) > max_distance:
                            mols.append((last_end, ci - 1))
                            last_end = ci
                    mols.append((last_end, len(cache) - 1))
                    for m in range(len(mols)):
                        starts = []
                        sizes = []
                        reads = []
                        bx = {}
                        bx_name = {}
                        mqual = 0
                        mori = 0
                        mchr = {}
                        for ix in range(mols[m][0], mols[m][1] + 1):
                            if bin(cache[ix][1])[2:].zfill(12)[::-1][10] == "1":
                                continue
                            mchr[cache[ix][2]] = mchr.get(cache[ix][2], 0) + 1
                            cigar = re.findall(r"(\d+[MIDNSHPX=])", cache[ix][5])
                            a_len = 0
                            for v in cigar:
                                if re.match(r"(\d+)[MDNX=]", v):
                                    a_len += int(re.match(r"(\d+)[MDNX=]", v).group(1))
                            if not starts:
                                starts.append(cache[ix][3] - 1
