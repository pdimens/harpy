import os
import re

localrules: all, concat_results
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

VERSION     = config.get("Workflow", {}).get('harpy-version', 'latest')
lr_platform = config.get("Workflow", {}).get("linkedreads", {}).get("type", 'none')
fqlist      = config["Inputs"]
bn_r        = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
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
        temp("{sample}.F.log")
    params:
        lr_platform
    shell: 
        "harpy-utils check-fastq {params} {input} > {output}"

rule check_reverse:
    input:
        get_fq2
    output:
        temp("{sample}.R.log")
    params:
        lr_platform
    shell: 
        "harpy-utils check-fastq {params} {input} > {output}"

rule concat_results:
    input:
        collect("{sample}.{FR}.log", sample = samplenames, FR = ["F","R"])
    output:
        "validate.fastq.tsv"
    shell:
        """
        {{
            echo -e "file\treads\tnoBX\tbadBX\tbadSamSpec\tbxNotLast"
            cat {input} | sort -k1
        }} > {output}
        """

rule create_report:
    input:
        data = "validate.fastq.tsv",
        ipynb = "workflow/validate_fastq.ipynb"
    output:
        tmp = temp("validate.fastq.tmp.ipynb"),
        ipynb = "validate.fastq.ipynb"
    params:
        lr_platform - lr_platform,
        infile = "-p infile " + os.path.abspath("validate.fastq.tsv")
    log:
        "logs/report.log"
    shell:
        """
        {{
            papermill -k python3 --cwd . --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params.infile}
            process-notebook {output.tmp} {params.lr_platform}
        }} 2> {log} > {output.ipynb}
        """

rule all:
    default_target: True
    input:
        "validate.fastq.ipynb"

