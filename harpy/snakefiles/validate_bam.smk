import os
import re
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["Workflow"]["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

lr_platform = config["Workflow"]["linkedreads"]["type"]
bamlist = config["Inputs"]
bamdict = dict(zip(bamlist, bamlist))
samplenames = {Path(i).stem for i in bamlist}

def get_alignments(wildcards):
    """returns a list with the bam file for the sample based on wildcards.sample"""
    r = re.compile(fr".*/({wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0]

rule check_bam:
    input:
        get_alignments
    output:
        temp("{sample}.log")
    params:
        lr_platform
    shell: 
        "check_bam {params} {input} > {output}"

rule concat_results:
    input:
        collect("{sample}.log", sample = samplenames)
    output:
        "validate.bam.tsv"
    shell:
        """
        {{
        echo -e "file\talignments\tnameMismatch\tnoMI\tnoBX\tbxNotLast\tbadBX"
        cat {input} | sort -k1
        }} > {output}
        """

rule create_report:
    input:
        data = "validate.bam.tsv",
        ipynb = "workflow/report/validate_bam.ipynb"
    output:
        tmp = temp("validate.bam.tmp.ipynb"),
        ipynb = "validate.bam.ipynb"
    params:
        lr_platform = lr_platform,
        static = "--no-progress-bar --log-level ERROR -k ir",
        sed_replace = 's/"injected-parameters"/"injected-parameters",\\n"remove-cell"/g'
    log:
        "logs/report.log"
    conda:
        "envs/report.yaml"
    container:
        "docker://pdimens/harpy:report_latest"
    shell:
        """
        {{
            papermill {params.static} {input.ipynb} {output.tmp} -p infile {input.data} -p platform {params.lr_platform}
            sed '{params.sed_replace}' {output.tmp}
        }} 2> {log} > {output.ipynb}
        """

rule all:
    default_target: True
    input:
        "validate.bam.ipynb"
