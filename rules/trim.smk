import os
import re
import glob

maxlen 		= config["maxlen"]
extra 		= config.get("extra", "") 
seq_dir 	= config["seq_directory"]
#fqext 		= config["fqext"]
#samplenames = config["samplenames"]

#flist = os.listdir(seq_dir)
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


rule trimFastp:
    input:
        fw   = get_fq1,
        rv   = get_fq2
    output:
        fw   = "Trim/{sample}.R1.fq.gz",
        rv   = "Trim/{sample}.R2.fq.gz",
        json = "Trim/logs/json/{sample}.fastp.json"
    log:
        html = "Trim/reports/{sample}.html",
        serr = "Trim/logs/err/{sample}.log"
    benchmark:
        "Benchmark/Trim/{sample}.txt"
    message:
        "Removing adapters + quality trimming: {wildcards.sample}"
    #wildcard_constraints: 
    #    sample = "[a-zA-Z0-9_-.]*"
    threads: 2
    params:
        maxlen = f"--max_len1 {maxlen}",
        extra = extra
    shell: 
        "fastp --trim_poly_g --cut_right --detect_adapter_for_pe {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {log.html} -j {output.json} 2> {log.serr}"

rule count_beadtags:
    input:
        "Trim/{sample}.R1.fq.gz"
    output: 
        temp("Trim/bxcount/{sample}.count.log")
    #wildcard_constraints:
    #    sample = "[a-zA-Z0-9_-.]*"
    message:
        "Counting barcode frequency: {wildcards.sample}"
    shell:
        "countBX.py {input} > {output}"

rule beadtag_counts_summary:
    input: 
        countlog = expand("Trim/bxcount/{sample}.count.log", sample = samplenames)
    output:
        "Trim/summary.bx.valid.html"
    message:
        "Summarizing sample barcode validation"
    script:
        "reportBxCount.Rmd"


rule log_runtime:
    output:
        "Trim/logs/harpy.trim.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    params:
        maxlen = f"--max_len1 {maxlen}",
        extra = extra
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy trim module ran using these parameters:\n\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n")
            _ = f.write("fastp trimming ran using:\n")
            _ = f.write("\tfastp --trim_poly_g --cut_right --detect_adapter_for_pe" + " ".join(params) + "\n")

rule createReport:
    input: 
        expand("Trim/logs/json/{sample}.fastp.json", sample = samplenames),
        expand("Trim/{sample}.{FR}.fq.gz", FR = ["R1", "R2"], sample = samplenames),
        "Trim/summary.bx.valid.html",
        "Trim/logs/harpy.trim.log"
    output:
        "Trim/trim.report.html"
    message:
        "Sequencing quality filtering and trimming is complete!"
    default_target: True
    shell: 
        "multiqc Trim/logs/json -m fastp --force --filename {output} --quiet --no-data-dir 2>/dev/null"