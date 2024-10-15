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

fqlist = config["inputs"]
outdir = config["output_directory"]
envdir      = os.getcwd() + "/.harpy_envs"
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}

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
    log:
        logfile = outdir + "/logs/report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/preflight_fastq.Rmd"

rule workflow_summary:
    default_target: True
    input:
        outdir + "/filecheck.fastq.html"
    run:
        os.makedirs(f"{outdir}/workflow/", exist_ok= True)
        summary_template = f"""
The harpy preflight fastq workflow ran using these parameters:

Validations were performed with:
    check_fastq.py sample.fastq > sample.txt

The Snakemake workflow was called via command line:
    {config["workflow_call"]}
"""
        with open(outdir + "/workflow/preflight.fastq.summary", "w") as f:
            f.write(summary_template)
