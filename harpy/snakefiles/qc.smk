containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)
wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

fqlist       = config["inputs"]
outdir       = config["output_directory"]
envdir       = os.path.join(os.getcwd(), outdir, "workflow", "envs")
min_len 	 = config["min_len"]
max_len 	 = config["max_len"]
extra 	     = config.get("extra", "") 
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
            fw   = outdir + "/{sample}.R1.fq.gz",
            rv   = outdir + "/{sample}.R2.fq.gz",
            html = outdir + "/reports/{sample}.html",
            json = outdir + "/reports/data/fastp/{sample}.fastp.json"
        log:
            serr = outdir + "/logs/fastp/{sample}.log"
        params:
            minlen = f"--length_required {min_len}",
            maxlen = f"--max_len1 {max_len}",
            trim_adapters = trim_arg,
            dedup = "-D" if dedup else "",
            extra = extra
        threads:
            workflow.cores
        conda:
            f"{envdir}/qc.yaml"
        shell: 
            """
            fastp --trim_poly_g --cut_right {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {output.html} -j {output.json} -R "{wildcards.sample} QC Report" 2> {log.serr}
            """
else:
    rule fastp:
        priority: 100
        input:
            fw   = get_fq1,
            rv   = get_fq2
        output:
            fq   = temp(outdir + "/fastp/{sample}.fastq"),
            html = outdir + "/reports/{sample}.html",
            json = outdir + "/reports/data/fastp/{sample}.fastp.json"
        log:
            serr = outdir + "/logs/fastp/{sample}.log"
        params:
            minlen = f"--length_required {min_len}",
            maxlen = f"--max_len1 {max_len}",
            trim_adapters = trim_arg,
            dedup = "-D" if dedup else "",
            extra = extra
        threads:
            workflow.cores
        conda:
            f"{envdir}/qc.yaml"
        shell: 
            """
            fastp --trim_poly_g --cut_right {params} --thread {threads} -i {input.fw} -I {input.rv} --stdout -h {output.html} -j {output.json} -R "{wildcards.sample} QC Report" 2> {log.serr} > {output.fq}
            """

    rule deconvolve:
        input:
            outdir + "/fastp/{sample}.fastq"
        output:
            temp(outdir + "/{sample}.fastq")
        log:
            outdir + "/logs/deconvolve/{sample}.deconvolve.log"
        params:
            kmer    = f"-k {decon_k}",
            windows = f"-w {decon_w}",
            density = f"-d {decon_d}",
            dropout = f"-a {decon_a}"
        threads:
            workflow.cores
        conda:
            f"{envdir}/qc.yaml"
        shell:
            "QuickDeconvolution -t {threads} -i {input} -o {output} {params} > {log} 2>&1"

    rule extract_forward:
        input:
            outdir + "/{sample}.fastq"
        output:
            outdir + "/{sample}.R1.fq.gz"
        params:
            "-1"
        container:
            None
        shell:
            "seqtk seq {params} {input} | gzip > {output}"

    use rule extract_forward as extract_reverse with:
        output:
            outdir + "/{sample}.R2.fq.gz"
        params:
            "-2"

rule check_barcodes:
    input:
        outdir + "/{sample}.R1.fq.gz"
    output: 
        temp(outdir + "/logs/bxcount/{sample}.count.log")
    container:
        None
    shell:
        "count_bx.py {input} > {output}"

rule barcode_report:
    input: 
        data = collect(outdir + "/logs/bxcount/{sample}.count.log", sample = samplenames),
        qmd = f"{outdir}/workflow/report/bx_count.qmd"
    output:
        report = f"{outdir}/reports/barcode.summary.html",
        qmd = temp(f"{outdir}/reports/barcode.summary.qmd")
    params:
        f"{outdir}/logs/bxcount/"
    log:
        outdir + "/logs/barcode.report.log"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp {input.qmd} {output.qmd}
        INPATH=$(realpath {params})
        quarto render {output.qmd} -l {log} --quiet -P indir:$INPATH
        """
   
rule qc_report:
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
    shell: 
        "multiqc {params} --filename {output}"

rule workflow_summary:
    default_target: True
    input:
        fq = collect(outdir + "/{sample}.{FR}.fq.gz", FR = ["R1", "R2"], sample = samplenames),
        bx_report = outdir + "/reports/barcode.summary.html" if not skip_reports else [],
        agg_report = outdir + "/reports/qc.report.html" if not skip_reports else []    
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
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/qc.summary", "w") as f:
            f.write("\n\n".join(summary))
