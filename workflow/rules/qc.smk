import os
import re
import sys
import glob
from rich.panel import Panel
from rich import print as rprint

maxlen 	  = config["maxlen"]
extra 	  = config.get("extra", "") 
seq_dir   = config["seq_directory"]
adapters  = config["adapters"]
skipreports = config["skipreports"]

flist = [os.path.basename(i) for i in glob.iglob(f"{seq_dir}/*") if not os.path.isdir(i)]
r = re.compile(r".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
fqlist = list(filter(r.match, flist))
bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onsuccess:
    print("")
    rprint(
        Panel(
            "The workflow has finished successfully! Find the results in [bold]QC/[/bold]",
            title = "[bold]harpy qc",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy qc",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )



def get_fq1(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    r = re.compile(r".*[\_\.][FR][1]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, lst))
    return fqlist

def get_fq2(wildcards):
    # code that returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    r = re.compile(r".*[\_\.][R][2]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
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
    params:
        maxlen = f"--max_len1 {maxlen}",
        tim_adapters = "--detect_adapter_for_pe" if adapters else "--disable_adapter_trimming",
        extra = extra
    threads:
        2
    conda:
        os.getcwd() + "/harpyenvs/qc.yaml"
    message:
        "Removing adapters + quality trimming: {wildcards.sample}"
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
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    message:
        "Summarizing sample barcode validation"
    script:
        "report/BxCount.Rmd"

rule log_runtime:
    output:
        "QC/workflow/qc.workflow.summary"
    params:
        maxlen = f"--max_len1 {maxlen}",
        tim_adapters = "--detect_adapter_for_pe" if adapters else "--disable_adapter_trimming",
        extra = extra
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy qc module ran using these parameters:\n\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n")
            _ = f.write("fastp trimming ran using:\n")
            _ = f.write("    fastp --trim_poly_g --cut_right " + " ".join(params) + "\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")

results = list()
results.append(expand("QC/logs/json/{sample}.fastp.json", sample = samplenames))
results.append(expand("QC/{sample}.{FR}.fq.gz", FR = ["R1", "R2"], sample = samplenames))
results.append("QC/workflow/qc.workflow.summary")
if not skipreports:
    results.append("QC/logs/barcode.summary.html")
    
rule createReport:
    default_target: True
    input: 
        results
    output:
        "QC/logs/qc.report.html"
    conda:
        os.getcwd() + "/harpyenvs/qc.yaml"
    message:
        "Sequencing quality filtering and trimming is complete!"
    shell: 
        """
        multiqc QC/logs/json -m fastp --force --filename {output} --quiet --title "QC Summary" --comment "This report aggregates trimming and quality control metrics reported by fastp" --no-data-dir 2>/dev/null
        """