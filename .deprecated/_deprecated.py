
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
        target = Path(seqfile).resolve()
        linkedfile = f"{linkdir}/{bn}.{FR}.fq{gz_ext}"
        try:
            _ = Path(linkedfile).symlink_to(target)
        except:
            pass
    return samplenames