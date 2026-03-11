import os
import re
from pathlib import Path

localrules: all, concat_results
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

VERSION     = config.get("Workflow", {}).get('harpy-version', 'latest')
lr_platform = config.get("Workflow", {}).get("linkedreads", {}).get("type", 'none')
bamlist     = config["Inputs"]
bamdict     = dict(zip(bamlist, bamlist))
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
        "harpy-utils check-bam {params} {input} > {output}"

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
        ipynb = "workflow/validate_bam.ipynb"
    output:
        tmp = temp("validate.bam.tmp.ipynb"),
        ipynb = "validate.bam.ipynb"
    params:
        lr_platform = lr_platform,
        infile = "-p infile " + os.path.abspath("validate.bam.tsv")
    log:
        "logs/report.log"
    shell:
        """
        {{
            papermill -k python3 --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params.infile}
            process-notebook {output.tmp} {params.lr_platform}
        }} 2> {log} > {output.ipynb}
        """

rule all:
    default_target: True
    input:
        "validate.bam.ipynb"
