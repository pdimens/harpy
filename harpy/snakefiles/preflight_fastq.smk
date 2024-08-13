containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import logging as pylogging
from pathlib import Path

fqlist = config["inputs"]
outdir = config["output_directory"]
envdir      = os.getcwd() + "/.harpy_envs"
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onstart:
    extra_logfile_handler = pylogging.FileHandler(snakemake_log)
    logger.logger.addHandler(extra_logfile_handler)

def get_fq1(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    #samples_FR = [i for i in fqlist if wildcards.sample in i]
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]1|[_\.]F|[_\.]R1(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

def get_fq2(wildcards):
    # returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    #samples_FR = [i for i in fqlist if wildcards.sample in i]
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]2|[_\.]R|[_\.]R2(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

rule check_forward:
    input:
        get_fq1
    output:
        temp(outdir + "/{sample}.F.log")
    container:
        None
    shell: 
        "check_fastq.py {input} > {output}"

rule check_reverse:
    input:
        get_fq2
    output:
        temp(outdir + "/{sample}.R.log")
    container:
        None
    shell: 
        "check_fastq.py {input} > {output}"

rule concat_results:
    input:
        collect(outdir + "/{sample}.{FR}.log", sample = samplenames, FR = ["F","R"])
    output:
        tmp = temp(outdir + "/filecheck.tmp"),
        final = outdir + "/filecheck.fastq.tsv"
    container:
        None
    shell:
        """
        cat {input} | sort -k1 > {output.tmp}
        echo -e "file\treads\tnoBX\tbadBX\tbadSamSpec\tbxNotLast\n$(cat {output.tmp})" > {output.final}
        """

rule create_report:
    input:
        outdir + "/filecheck.fastq.tsv"
    output:
        outdir + "/filecheck.fastq.html"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/preflight_fastq.Rmd"

rule workflow_summary:
    default_target: True
    input:
        outdir + "/filecheck.fastq.html"
    message:
        "Summarizing the workflow: {output}"
    run:
        os.makedirs(f"{outdir}/workflow/", exist_ok= True)
        with open(outdir + "/workflow/preflight.fastq.summary", "w") as f:
            _ = f.write("The harpy preflight fastq workflow ran using these parameters:\n\n")
            _ = f.write("validations were performed with:\n")
            _ = f.write("    check_fastq.py sample.fastq > sample.txt\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")