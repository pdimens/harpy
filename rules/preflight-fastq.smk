import os
import re
import glob

seq_dir = config["seq_directory"]
out_dir = f"{seq_dir}/Preflight/"

flist = [os.path.basename(i) for i in glob.iglob(f"{seq_dir}/*") if not os.path.isdir(i)]
r = re.compile(".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
fqlist = list(filter(r.match, flist))
bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])

def get_fq1(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    r = re.compile(".*[\_\.][FR][1]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, lst))
    return fqlist

def get_fq2(wildcards):
    # code that returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    r = re.compile(".*[\_\.][R][2]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, lst))
    return fqlist

rule checkForward:
    input:
        get_fq1
    output:
        temp(out_dir + "{sample}.F.log")
    message:
        "Processing forward reads: {wildcards.sample}"
    shell: 
        "checkFASTQ.py {input} > {output}"

rule checkReverse:
    input:
        get_fq2
    output:
        temp(out_dir + "{sample}.R.log")
    message:
        "Processing reverse reads: {wildcards.sample}"
    shell: 
        "checkFASTQ.py {input} > {output}"

rule mergeChecks:
    default_target: True
    input:
        expand(out_dir + "{sample}.{FR}.log", sample = samplenames, FR = ["F","R"])
    output:
        tmp = temp(out_dir + "validations.tmp"),
        final = out_dir + "validations.fastq.tsv"
    message:
        "Concatenating results"
    shell:
        """
        cat {input} | sort -k1 > {output.tmp}
        echo -e "file\treads\tnoBX\tbadBX\tbadSamSpec\n$(cat {output.tmp})" > {output.final}
        """


# TODO add rmd report
#rule createReport:

