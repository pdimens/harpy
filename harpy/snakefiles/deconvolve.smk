containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import logging as pylogging

envdir      = os.getcwd() + "/.harpy_envs"
fqlist      = config["inputs"]
outdir      = config["output_directory"]
kmer_length = config["kmer_length"]
window_size = config["window_size"]
density 	= config["density"] 
dropout     = config["dropout"]
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
        with open(outdir + "/workflow/deconvolve.summary", "w") as f:
            _ = f.write("The harpy deconvolve workflow ran using these parameters:\n\n")
            _ = f.write("fastq files were interleaved with seqtk:\n")
            _ = f.write("    seqtk mergepe forward.fq reverse.fq\n")
            _ = f.write("Deconvolution occurred using QuickDeconvolution:\n")
            _ = f.write(f"   QuickDeconvolution -t threads -i infile.fq -o output.fq -k {kmer_length} -w {window_size} -d {density} -a {dropout}\n")
            _ = f.write("The interleaved output was split back into forward and reverse reads with seqtk:\n")
            _ = f.write("    seqtk -1 interleaved.fq | gzip > file.R1.fq.gz\n")
            _ = f.write("    seqtk -2 interleaved.fq | gzip > file.R2.fq.gz\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")
