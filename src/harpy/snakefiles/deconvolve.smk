containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import logging as pylogging
from datetime import datetime
from rich.panel import Panel
from rich import print as rprint

envdir      = os.getcwd() + "/.harpy_envs"
fqlist      = config["inputs"]
outdir      = config["output_directory"]
kmer_length = config["kmer_length"]
window_size = config["window_size"]
density 	= config["density"] 
dropout     = config["dropout"]

bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onstart:
    os.makedirs(f"{outdir}/logs/snakemake", exist_ok = True)
    dt_string = datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
    extra_logfile_handler = pylogging.FileHandler(f"{outdir}/logs/snakemake/{dt_string}.snakelog")
    logger.logger.addHandler(extra_logfile_handler)

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]",
            title = "[bold]harpy deconvolve",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file for more details:\n[bold]{outdir}/logs/snakemake/{dt_string}.snakelog[/bold]",
            title = "[bold]harpy deconvolve",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

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
    message:
        "Interleaving: {wildcards.sample}"
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

rule workflow_summary:
    default_target: True
    input:
        collect(outdir + "/{sample}.{FR}.fq.gz", FR = ["R1", "R2"], sample = samplenames),
    message:
        "Summarizing the workflow: {output}"
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
