import os
import re

wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

lr_type       = config["Workflow"]["linkedreads"]["type"]
skip_reports  = config["Workflow"]["reports"]["skip"]
fqlist        = config["Inputs"]
min_len 	  = config["Parameters"]["min-len"]
max_len 	  = config["Parameters"]["max-len"]
extra 	      = config["Parameters"].get("extra", "") 
trim_adapters = config["Parameters"].get("trim_adapters", None)
dedup         = config["Parameters"]["deduplicate"]
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
if trim_adapters:
    trim_arg = "--detect_adapter_for_pe" if trim_adapters == "auto" else f"--adapter_fasta {trim_adapters}"
else:
    trim_arg = "--disable_adapter_trimming"

def get_fq1(wildcards):
    '''returns a list of fastq files for read 1 based on *wildcards.sample*'''
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]1|[_\.]F|[_\.]R1(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

def get_fq2(wildcards):
    '''returns a list of fastq files for read 2 based on *wildcards.sample*'''
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]2|[_\.]R|[_\.]R2(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

rule fastp:
    priority: 100
    input:
        fw   = get_fq1,
        rv   = get_fq2
    output:
        fw   = "{sample}.R1.fq.gz",
        rv   = "{sample}.R2.fq.gz",
        html = "reports/{sample}.html",
        json = "reports/data/fastp/{sample}.fastp.json"
    log:
        serr = "logs/fastp/{sample}.log"
    params:
        static = "--trim_poly_g --cut_right",
        minlen = f"--length_required {min_len}",
        maxlen = f"--max_len1 {max_len}",
        trim_adapters = trim_arg,
        dedup = "-D" if dedup else "--dont_eval_duplication",
        title = lambda wc: f"-R \"{wc.sample} QC Report\"",
        extra = extra
    threads:
        workflow.cores
    conda:
        "envs/qc.yaml"
    container:
        "docker://pdimens/harpy:qc_dev"
    shell: 
        "fastp {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {output.html} -j {output.json} 2> {log.serr}"

rule barcode_stats:
    input:
        "{sample}.R1.fq.gz"
    output: 
        temp("logs/bxcount/{sample}.count.log")
    params:
        lr_type
    shell:
        "count_bx {params} {input} > {output}"

rule barcode_report:
    input:
        data = collect("logs/bxcount/{sample}.count.log", sample = samplenames),
        ipynb = f"workflow/qc_bx_stats.ipynb"
    output:
        tmp = temp("reports/barcode.summary.tmp.ipynb"),
        ipynb = "reports/barcode.summary.ipynb"
    log:
        "logs/barcode.report.log"
    params:
        indir = "-p indir " + os.path.abspath("logs/bxcount"),
        lr = lr_type
    shell:
        """
        {{
            papermill --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params.indir}
            process_notebook {params.lr} {output.tmp}
        }} 2> {log} > {output.ipynb}
        """

rule qc_report:
    input: 
        collect("reports/data/fastp/{sample}.fastp.json", sample = samplenames)
    output:
        "reports/qc.report.html"
    log:
        "logs/multiqc.log"
    params:
        module = "-m fastp",
        options = "-n stdout --no-ai --no-version-check --force --quiet --no-data-dir",
        title = "--title \"QC Summary\"",
        comment = "--comment \"This report aggregates trimming and quality control metrics reported by fastp.\"",
        logdir = "reports/data/fastp/"
    conda:
        "envs/qc.yaml"
    container:
        "docker://pdimens/harpy:qc_dev"
    shell: 
        "multiqc {params} > {output} 2> {log}"

rule all:
    default_target: True
    input:
        fq = collect("{sample}.{FR}.fq.gz", FR = ["R1", "R2"], sample = samplenames),
        bx_report = "reports/barcode.summary.ipynb" if not skip_reports and lr_type != "none" else [],
        agg_report = "reports/qc.report.html" if not skip_reports else []    
