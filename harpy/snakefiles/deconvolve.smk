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
    sample = r"[a-zA-Z0-9._-]+"

fqlist      = config["inputs"]
outdir      = config["output_directory"]
envdir      = os.path.join(os.getcwd(), outdir, "workflow", "envs")
kmer_length = config["kmer_length"]
window_size = config["window_size"]
density 	= config["density"] 
dropout     = config["dropout"]
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}

def get_fq1(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]1|[_\.]F|[_\.]R1(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

def get_fq2(wildcards):
    # returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]2|[_\.]R|[_\.]R2(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

rule interleave:
    input:
        fw = get_fq1,
        rv = get_fq2
    output: 
        temp(outdir + "/workflow/input/{sample}.fastq"),
    container:
        None
    shell:
        "seqtk mergepe {input} > {output}"

rule deconvolve:
    input:
        outdir + "/workflow/input/{sample}.fastq"
    output:
        temp(outdir + "/{sample}.fastq")
    log:
        outdir + "/logs/{sample}.deconvolve.log"
    params:
        kmer    = f"-k {kmer_length}",
        windows = f"-w {window_size}",
        density = f"-d {density}",
        dropout = f"-a {dropout}"
    threads:
        2
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

rule workflow_summary:
    default_target: True
    input:
        collect(outdir + "/{sample}.{FR}.fq.gz", FR = ["R1", "R2"], sample = samplenames),
    run:
        summary = ["The harpy deconvolve workflow ran using these parameters:"]
        interleave = "fastq files were interleaved with seqtk:\n"
        interleave += "\tseqtk mergepe forward.fq reverse.fq"
        summary.append(interleave)
        deconv = "Deconvolution occurred using QuickDeconvolution:\n"
        deconv += f"\tQuickDeconvolution -t threads -i infile.fq -o output.fq -k {kmer_length} -w {window_size} -d {density} -a {dropout}"
        summary.append(deconv)
        recover = "The interleaved output was split back into forward and reverse reads with seqtk:\n"
        recover += "\tseqtk -1 interleaved.fq | gzip > file.R1.fq.gz\n"
        recover += "\tseqtk -2 interleaved.fq | gzip > file.R2.fq.gz"
        summary.append(recover)
        sm = "Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/deconvolve.summary", "w") as f:  
            f.write("\n\n".join(summary))