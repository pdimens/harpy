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
ignore_bx    = config["ignore_bx"]
trim_adapters = config.get("trim_adapters", None)
dedup        = config["deduplicate"]
deconvolve   = config.get("deconvolve", False)
if deconvolve:
    decon_k = deconvolve["kmer_length"]
    decon_w = deconvolve["window_size"]
    decon_d = deconvolve["density"]
    decon_a = deconvolve["dropout"]
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

if not deconvolve:
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
else:
    rule fastp:
        priority: 100
        input:
            fw   = get_fq1,
            rv   = get_fq2
        output:
            fq   = temp("fastp/{sample}.fastq"),
            html = "reports/{sample}.html",
            json = "reports/data/fastp/{sample}.fastp.json"
        log:
            serr = "logs/fastp/{sample}.log"
        params:
            static = "--trim_poly_g --cut_right",
            minlen = f"--length_required {min_len}",
            maxlen = f"--max_len1 {max_len}",
            trim_adapters = trim_arg,
            dedup = "-D" if dedup else "",
            title = lambda wc: f"-R \"{wc.sample} QC Report\"",
            extra = extra
        threads:
            workflow.cores
        conda:
            "envs/qc.yaml"
        shell: 
            "fastp {params} --thread {threads} -i {input.fw} -I {input.rv} --stdout -h {output.html} -j {output.json} 2> {log.serr} > {output.fq}"

    rule deconvolve:
        input:
            "fastp/{sample}.fastq"
        output:
            temp("{sample}.fastq")
        log:
            "logs/deconvolve/{sample}.deconvolve.log"
        params:
            kmer    = f"-k {decon_k}",
            windows = f"-w {decon_w}",
            density = f"-d {decon_d}",
            dropout = f"-a {decon_a}"
        threads:
            workflow.cores
        conda:
            "envs/qc.yaml"
        shell:
            "QuickDeconvolution -t {threads} -i {input} -o {output} {params} > {log} 2>&1"

    rule extract_forward:
        input:
            "{sample}.fastq"
        output:
            "{sample}.R1.fq.gz"
        params:
            "-1"
        container:
            None
        shell:
            "seqtk seq {params} {input} | gzip > {output}"

    use rule extract_forward as extract_reverse with:
        output:
            "{sample}.R2.fq.gz"
        params:
            "-2"

rule barcode_stats:
    input:
        "{sample}.R1.fq.gz"
    output: 
        temp("logs/bxcount/{sample}.count.log")
    container:
        None
    shell:
        "count_bx {input} > {output}"

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
        qmd = f"workflow/report/bx_count.qmd"
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
        quarto render {output.qmd} --log {log} --quiet -P indir:$INPATH
        """
   
rule qc_report:
    input: 
        collect("reports/data/fastp/{sample}.fastp.json", sample = samplenames)
    output:
        "reports/qc.report.html"
    log:
        "logs/multiqc.log"
    params:
        logdir = "reports/data/fastp/",
        module = "-m fastp",
        options = "--no-ai --no-version-check --force --quiet --no-data-dir",
        title = "--title \"QC Summary\"",
        comment = "--comment \"This report aggregates trimming and quality control metrics reported by fastp.\""
    conda:
        "envs/qc.yaml"
    shell: 
        "multiqc {params} --filename {output} 2> {log}"

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
        if deconvolve:
            deconv = "Deconvolution occurred using QuickDeconvolution:\n"
            deconv += "\tQuickDeconvolution -t threads -i infile.fq -o output.fq -k {decon_k} -w {decon_w} -d {decon_d} -a {decon_a}"
            summary.append(deconv)
            interlv = "The interleaved output was split back into forward and reverse reads with seqtk:\n"
            interlv += "\tseqtk -1 interleaved.fq | gzip > file.R1.fq.gz\n"
            interlv += "\tseqtk -2 interleaved.fq | gzip > file.R2.fq.gz"
            summary.append(interlv)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/qc.summary", "w") as f:
            f.write("\n\n".join(summary))
