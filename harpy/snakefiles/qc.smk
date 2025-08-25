containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

fqlist       = config["inputs"]
min_len 	 = config["min_len"]
max_len 	 = config["max_len"]
extra 	     = config.get("extra", "") 
lr_type    = config["linkedread_type"]
ignore_bx = lr_type == "none"
trim_adapters = config.get("trim_adapters", None)
dedup        = config["deduplicate"]
skip_reports  = config["reports"]["skip"]
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
if trim_adapters:
    trim_arg = "--detect_adapter_for_pe" if trim_adapters == "auto" else f"--adapter_fasta {trim_adapters}"
else:
    trim_arg = "--disable_adapter_trimming"

def get_fq1(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]1|[_\.]F|[_\.]R1(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

def get_fq2(wildcards):
    # returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
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
    shell: 
        "fastp {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {output.html} -j {output.json} 2> {log.serr}"

rule barcode_stats:
    input:
        "{sample}.R1.fq.gz"
    output: 
        temp("logs/bxcount/{sample}.count.log")
    params:
        lr_type
    container:
        None
    shell:
        "count_bx {params} {input} > {output}"

rule configure_report:
    input:
        yaml = f"workflow/report/_quarto.yml",
        scss = f"workflow/report/_harpy.scss"
    output:
        yaml = temp("reports/_quarto.yml"),
        scss = temp("reports/_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule barcode_report:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        data = collect("logs/bxcount/{sample}.count.log", sample = samplenames),
        qmd = f"workflow/report/qc_bx_stats.qmd"
    output:
        report = "reports/barcode.summary.html",
        qmd = temp("reports/barcode.summary.qmd")
    params:
        f"logs/bxcount/"
    log:
        "logs/barcode.report.log"
    conda:
        "envs/r.yaml"
    retries:
        3
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        INPATH=$(realpath {params})
        quarto render {output.qmd} --no-cache --log {log} --quiet -P indir:$INPATH
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
    shell: 
        "multiqc {params} > {output} 2> {log}"

rule workflow_summary:
    default_target: True
    input:
        fq = collect("{sample}.{FR}.fq.gz", FR = ["R1", "R2"], sample = samplenames),
        bx_report = "reports/barcode.summary.html" if not skip_reports and not ignore_bx else [],
        agg_report = "reports/qc.report.html" if not skip_reports else []    
    params:
        minlen = f"--length_required {min_len}",
        maxlen = f"--max_len1 {max_len}",
        trim_adapters = trim_arg,
        dedup = "-D" if dedup else "",
        extra = extra
    run:
        summary = ["The harpy qc workflow ran using these parameters:"]
        fastp = "fastp ran using:\n"
        fastp += f"\tfastp --trim_poly_g --cut_right {params}"
        summary.append(fastp)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/qc.summary", "w") as f:
            f.write("\n\n".join(summary))
