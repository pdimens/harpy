containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import glob
from pathlib import Path
from rich import print as rprint
from rich.panel import Panel

fqlist = config["inputs"]
out_dir = config["output_directory"]
envdir      = os.getcwd() + "/.harpy_envs"
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy preflight fastq",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{out_dir}[/bold]",
            title = "[bold]harpy preflight fastq",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

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
        temp(out_dir + "/{sample}.F.log")
    container:
        None
    message:
        "Processing forward reads: {wildcards.sample}"
    shell: 
        "check_fastq.py {input} > {output}"

rule check_reverse:
    input:
        get_fq2
    output:
        temp(out_dir + "/{sample}.R.log")
    message:
        "Processing reverse reads: {wildcards.sample}"
    container:
        None
    shell: 
        "check_fastq.py {input} > {output}"

rule merge_checks:
    input:
        collect(out_dir + "/{sample}.{FR}.log", sample = samplenames, FR = ["F","R"])
    output:
        tmp = temp(out_dir + "/filecheck.tmp"),
        final = out_dir + "/filecheck.fastq.tsv"
    container:
        None
    message:
        "Concatenating results"
    shell:
        """
        cat {input} | sort -k1 > {output.tmp}
        echo -e "file\treads\tnoBX\tbadBX\tbadSamSpec\tbxNotLast\n$(cat {output.tmp})" > {output.final}
        """

rule create_report:
    input:
        out_dir + "/filecheck.fastq.tsv"
    output:
        out_dir + "/filecheck.fastq.html"
    conda:
        f"{envdir}/r.yaml"
    message:
        "Producing report"
    script:
        "report/preflight_fastq.Rmd"

rule workflow_summary:
    default_target: True
    input:
        out_dir + "/filecheck.fastq.html"
    message:
        "Summarizing the workflow: {output}"
    run:
        os.makedirs(f"{out_dir}/workflow/", exist_ok= True)
        with open(out_dir + "/workflow/preflight.fastq.summary", "w") as f:
            _ = f.write("The harpy preflight fastq workflow ran using these parameters:\n\n")
            _ = f.write("validations were performed with:\n")
            _ = f.write("    check_fastq.py sample.fastq > sample.txt\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")