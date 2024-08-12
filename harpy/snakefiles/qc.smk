containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import logging as pylogging

envdir       = os.getcwd() + "/.harpy_envs"
fqlist       = config["inputs"]
outdir       = config["output_directory"]
min_len 	 = config["min_len"]
max_len 	 = config["max_len"]
extra 	     = config.get("extra", "") 
trimadapters = config["trim_adapters"]
dedup        = config["deduplicate"]
deconvolve   = config.get("deconvolve", False)
if deconvolve:
    decon_k = deconvolve["kmer_length"]
    decon_w = deconvolve["window_size"]
    decon_d = deconvolve["density"]
    decon_a = deconvolve["dropout"]
skipreports  = config["skip_reports"]
snakemake_log = config["snakemake_log"]

bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onstart:
    extra_logfile_handler = pylogging.FileHandler(snakemake_log)
    logger.logger.addHandler(extra_logfile_handler)

def get_fq1(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]1|[_\.]F|[_\.]R1(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

def get_fq2(wildcards):
    # returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]2|[_\.]R|[_\.]R2(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

if not deconvolve:
    rule qc_fastp:
        input:
            fw   = get_fq1,
            rv   = get_fq2
        output:
            fw   = outdir + "/{sample}.R1.fq.gz",
            rv   = outdir + "/{sample}.R2.fq.gz",
            json = outdir + "/reports/data/fastp/{sample}.fastp.json"
        log:
            html = outdir + "/reports/{sample}.html",
            serr = outdir + "/logs/fastp/{sample}.log"
        params:
            minlen = f"--length_required {min_len}",
            maxlen = f"--max_len1 {max_len}",
            trim_adapters = "--detect_adapter_for_pe" if trimadapters else "--disable_adapter_trimming" ,
            dedup = "-D" if dedup else "",
            extra = extra
        threads:
            workflow.cores
        conda:
            f"{envdir}/qc.yaml"
        message:
            "Quality trimming" + (", removing adapters" if not trimadapters else "") + (", removing PCR duplicates" if dedup else "") + ": {wildcards.sample}"
        shell: 
            """
            fastp --trim_poly_g --cut_right {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {log.html} -j {output.json} -R "{wildcards.sample} QC Report" 2> {log.serr}
            """
else:
    rule qc_fastp:
        input:
            fw   = get_fq1,
            rv   = get_fq2
        output:
            fq   = temp(outdir + "/fastp/{sample}.fastq"),
            json = outdir + "/reports/data/fastp/{sample}.fastp.json"
        log:
            html = outdir + "/reports/{sample}.html",
            serr = outdir + "/logs/fastp/{sample}.log"
        params:
            minlen = f"--length_required {min_len}",
            maxlen = f"--max_len1 {max_len}",
            trim_adapters = "--detect_adapter_for_pe" if trimadapters else "--disable_adapter_trimming" ,
            dedup = "-D" if dedup else "",
            extra = extra
        threads:
            workflow.cores
        conda:
            f"{envdir}/qc.yaml"
        message:
            "Quality trimming" + (", removing adapters" if not trimadapters else "") + (", removing PCR duplicates" if dedup else "") + ": {wildcards.sample}"
        shell: 
            """
            fastp --trim_poly_g --cut_right {params} --thread {threads} -i {input.fw} -I {input.rv} --stdout -h {log.html} -j {output.json} -R "{wildcards.sample} QC Report" 2> {log.serr} > {output.fq}
            """

    rule deconvolve:
        input:
            outdir + "/fastp/{sample}.fastq"
        output:
            temp(outdir + "/{sample}.fastq")
        log:
            outdir + "/logs/{sample}.deconvolve.log"
        params:
            kmer    = f"-k {decon_k}",
            windows = f"-w {decon_w}",
            density = f"-d {decon_d}",
            dropout = f"-a {decon_a}"
        threads:
            workflow.cores
        conda:
            f"{envdir}/qc.yaml"
        message:
            "Performing deconvolution: {wildcards.sample}"
        shell:
            "QuickDeconvolution -t {threads} -i {input} -o {output} {params} > {log} 2>&1"

    rule recover_forward:
        input:
            outdir + "/{sample}.fastq"
        output:
            outdir + "/{sample}.R1.fq.gz"
        params:
            "-1"
        container:
            None
        message:
            "Extracting deconvolved forward reads: {wildcards.sample}"
        shell:
            "seqtk seq {params} {input} | gzip > {output}"

    use rule recover_forward as recover_reverse with:
        output:
            outdir + "/{sample}.R2.fq.gz"
        params:
            "-2"

rule count_beadtags:
    input:
        outdir + "/{sample}.R1.fq.gz"
    output: 
        temp(outdir + "/logs/bxcount/{sample}.count.log")
    message:
        "Counting barcode frequency: {wildcards.sample}"
    container:
        None
    shell:
        "count_bx.py {input} > {output}"

rule beadtag_counts_summary:
    input: 
        countlog = collect(outdir + "/logs/bxcount/{sample}.count.log", sample = samplenames)
    output:
        outdir + "/reports/barcode.summary.html"
    conda:
        f"{envdir}/r.yaml"
    message:
        "Summarizing sample barcode validation"
    script:
        "report/bx_count.Rmd"
   
rule create_report:
    input: 
        collect(outdir + "/reports/data/fastp/{sample}.fastp.json", sample = samplenames)
    output:
        outdir + "/reports/qc.report.html"
    params:
        logdir = f"{outdir}/reports/data/fastp/",
        module = "-m fastp",
        options = "--no-version-check --force --quiet --no-data-dir",
        title = "--title \"QC Summary\"",
        comment = "--comment \"This report aggregates trimming and quality control metrics reported by fastp.\""
    conda:
        f"{envdir}/qc.yaml"
    message:
        "Aggregating fastp reports"
    shell: 
        "multiqc {params} --filename {output}"

rule workflow_summary:
    default_target: True
    input:
        fq = collect(outdir + "/{sample}.{FR}.fq.gz", FR = ["R1", "R2"], sample = samplenames),
        bx_report = outdir + "/reports/barcode.summary.html" if not skipreports else [],
        agg_report = outdir + "/reports/qc.report.html" if not skipreports else []    
    params:
        minlen = f"--length_required {min_len}",
        maxlen = f"--max_len1 {max_len}",
        trim_adapters = "--detect_adapter_for_pe" if trimadapters else "--disable_adapter_trimming" ,
        dedup = "-D" if dedup else "",
        extra = extra
    message:
        "Summarizing the workflow: {output}"
    run:
        with open(outdir + "/workflow/qc.summary", "w") as f:
            _ = f.write("The harpy qc workflow ran using these parameters:\n\n")
            _ = f.write("fastp trimming ran using:\n")
            _ = f.write("    fastp --trim_poly_g --cut_right " + " ".join(params) + "\n")
            if deconvolve:
                _ = f.write("Deconvolution occurred using QuickDeconvolution:\n")
                _ = f.write(f"   QuickDeconvolution -t threads -i infile.fq -o output.fq -k {decon_k} -w {decon_w} -d {decon_d} -a {decon_a}\n")
                _ = f.write("The interleaved output was split back into forward and reverse reads with seqtk:\n")
                _ = f.write("    seqtk -1 interleaved.fq | gzip > file.R1.fq.gz\n")
                _ = f.write("    seqtk -2 interleaved.fq | gzip > file.R2.fq.gz\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")