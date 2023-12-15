import os
import re
import glob

maxlen 	  = config["maxlen"]
extra 	  = config.get("extra", "") 
seq_dir   = config["seq_directory"]
adapters  = config["adapters"]

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
        fw   = "QC/{sample}.R1.fq.gz",
        rv   = "QC/{sample}.R2.fq.gz",
        json = "QC/logs/json/{sample}.fastp.json"
    log:
        html = "QC/logs/fastp_reports/{sample}.html",
        serr = "QC/logs/fastp_logs/{sample}.log"
    benchmark:
        ".Benchmark/QC/{sample}.txt"
    message:
        "Removing adapters + quality trimming: {wildcards.sample}"
    threads:
        2
    params:
        maxlen = f"--max_len1 {maxlen}",
        tim_adapters = "--detect_adapter_for_pe" if adapters else "--disable_adapter_trimming",
        extra = extra
    shell: 
        """
        fastp --trim_poly_g --cut_right {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {log.html} -j {output.json} -R "{wildcards.sample} QC Report" 2> {log.serr}
        """

rule count_beadtags:
    input:
        "QC/{sample}.R1.fq.gz"
    output: 
        temp("QC/logs/bxcount/{sample}.count.log")
    message:
        "Counting barcode frequency: {wildcards.sample}"
    shell:
        "countBX.py {input} > {output}"

rule beadtag_counts_summary:
    input: 
        countlog = expand("QC/logs/bxcount/{sample}.count.log", sample = samplenames)
    output:
        "QC/logs/barcode.summary.html"
    message:
        "Summarizing sample barcode validation"
    script:
        "reportBxCount.Rmd"

rule log_runtime:
    output:
        "QC/logs/harpy.qc.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    params:
        maxlen = f"--max_len1 {maxlen}",
        tim_adapters = "--detect_adapter_for_pe" if adapters else "--disable_adapter_trimming",
        extra = extra
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy qc module ran using these parameters:\n\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n")
            _ = f.write("fastp trimming ran using:\n")
            _ = f.write("    fastp --trim_poly_g --cut_right " + " ".join(params) + "\n")

rule createReport:
    default_target: True
    input: 
        expand("QC/logs/json/{sample}.fastp.json", sample = samplenames),
        expand("QC/{sample}.{FR}.fq.gz", FR = ["R1", "R2"], sample = samplenames),
        "QC/logs/barcode.summary.html",
        "QC/logs/harpy.qc.log"
    output:
        "QC/logs/qc.report.html"
    message:
        "Sequencing quality filtering and trimming is complete!"
    shell: 
        """
        multiqc QC/logs/json -m fastp --force --filename {output} --quiet --title "QC Summary" --comment "This report aggregates trimming and quality control metrics reported by fastp" --no-data-dir 2>/dev/null
        """